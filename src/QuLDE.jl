"""
    Based on :  arxiv.org/abs/1807.04553

    * Uses Taylor expansion

    x' = Ax + b

    * A - input matrix.
    * b - input vector.
    * x - inital vector
    * t - time to be evaluated at

    * Consists of two parts: A is unitary - QuLDEUnitParam() & A is non unitary - QuLDEnonUnitParam()
    * quldecircuit() - generates qunatum cicuit

"""

struct QuLDEUnitParam{CPType, L, Q, VType <: GeneralMatrixBlock{1}, SType <: GeneralMatrixBlock{Q}}
    k::Int # Taylor expansion upto k
    t::L # time
    C_tilda::CPType
    D_tilda::CPType
    N::CPType
    V::VType
    VS1::SType # Q is input size of GMblock
    VS2::SType # Q
    WS1::SType # Q
    WS2::SType # Q

    function QuLDEUnitParam(k::Int,t::L,prob::QuLDEProblem) where L
        CPType = eltype(prob.u0)
        C(m) = norm(prob.u0)*(opnorm(prob.A)*t)^(m)/factorial(m)
        D(m) = norm(prob.b)*(opnorm(prob.A)*t)^(m-1)*t/factorial(m)
        C_tilda = 0
        D_tilda = 0
        for i in 1:k
            C_tilda = C_tilda + C(i)
            D_tilda = D_tilda + D(i)
        end
        C_tilda = C_tilda + C(0)
        N = sqrt(C_tilda + D_tilda)
        C_tilda = sqrt(C_tilda)
        D_tilda = sqrt(D_tilda)
        V = CPType[C_tilda/N D_tilda/N; D_tilda/N -1*C_tilda/N]
        V = matblock(V)

        VS1 = rand(CPType,k+1,k+1)
        VS2 = rand(CPType,k+1,k+1)
        for j in 0:k
            VS1[j+1,1] = sqrt(C(j))/C_tilda
            VS2[j+1,1] = sqrt(D(j+1))/D_tilda
        end
        VS2[k+1,1] = 0;
        VS1 = -1*qr(VS1).Q
        VS2 = -1*qr(VS2).Q
        WS1 = convert(typeof(VS1),VS1')
        WS2 = convert(typeof(VS1),VS2')
        VS1 = matblock(VS1)
        VS2 = matblock(VS2)
        WS1 = matblock(WS1)
        WS2 = matblock(WS2)
        n = log2i(k+1)
        new{CPType, L, n, typeof(V), typeof(VS1)}(k,t,C_tilda,D_tilda,N,V,VS1,VS2,WS1,WS2)
    end
end

struct QuLDEnonUnitParam{CPType ,L, Q, NL, Nbit, VType <: GeneralMatrixBlock{1}, SType <: GeneralMatrixBlock{Q}, LType <: GeneralMatrixBlock{NL}}
    k::Int # Taylor expansion upto k
    t::L # time
    l::Int # input size of LType block set to 2
    C_tilda::CPType #G1 as per the article
    D_tilda::CPType #G2 as per the article
    N::CPType
    V::VType
    VS1::SType # Q is input size of GMblock
    VS2::SType # Q
    WS1::SType # Q
    WS2::SType # Q
    VT::LType # NL is input size of GMblock
    WT::LType # NL
    F::Array{GeneralMatrixBlock{Nbit,Nbit,CPType,Array{CPType,2}},1} # Nbit is input size

    function QuLDEnonUnitParam(k::Int,t::L,l::Int,prob::QuLDEProblem) where L
        CPType = eltype(prob.u0)
        C(m) = norm(prob.u0)*((opnorm(prob.A)*t*2)^(m))/factorial(m) # alphas equal 1/2
        D(n) = norm(prob.b)*(opnorm(prob.A)*t*2)^(n-1)*t/factorial(n)
        nbit, = size(prob.u0)
        nbit = log2i(nbit)
        Mu = prob.A/opnorm(prob.A)
        B1 = convert(Array{CPType,2},1/2*(Mu + Mu'))
        B2 = convert(Array{CPType,2},-im/2*(Mu - Mu'))
        iden = Matrix{CPType}(I,size(B1))
        F = Array{GeneralMatrixBlock{nbit,nbit,CPType,Array{CPType,2}},1}(undef,4)
        F[1] = matblock(B1 + im*sqrt(iden - B1*B1))
        F[2] = matblock(B1 - im*sqrt(iden - B1*B1))
        F[3] = matblock(im*B2 - sqrt(iden - B2*B2))
        F[4] = matblock(im*B2 + sqrt(iden - B2*B2))
        #tested
        C_tilda = 0
        D_tilda = 0
        for i in 1:k
            C_tilda = C_tilda + C(i)
            D_tilda = D_tilda + D(i)
        end
        C_tilda = C_tilda + C(0)
        N = sqrt(C_tilda+ D_tilda)
        C_tilda = sqrt(C_tilda)
        D_tilda = sqrt(D_tilda)
        V = CPType[C_tilda/N D_tilda/N; D_tilda/N -1*C_tilda/N]
        V = matblock(V)
        #tested
        VS1 = rand(CPType,2^k,2^k)
        VS2 = rand(CPType,2^k,2^k)
        VS1[:,1] = zero(VS1[:,1])
        VS2[:,1] = zero(VS2[:,1])
        for j in 0:k
            VS1[(2^k - 2^(k-j) + 1),1] = sqrt(C(j))/C_tilda
            VS2[(2^k - 2^(k-j) + 1),1] = sqrt(D(j+1))/D_tilda
        end
        VS2[2^k,1] = 0;
        VS1 = -1*qr(VS1).Q
        VS2 = -1*qr(VS2).Q
        WS1 = convert(typeof(VS1),VS1')
        WS2 = convert(typeof(VS1),VS2')
        VS1 = matblock(VS1)
        VS2 = matblock(VS2)
        WS1 = matblock(WS1)
        WS2 = matblock(WS2)

        VT = rand(CPType,2^l,2^l)
        VT[:,1] = 0.5*ones(2^l,1)
        VT = -1*qr(VT).Q
        WT = convert(typeof(VT),VT')
        VT = matblock(VT)
        WT = matblock(WT)
        new{CPType, L, k, l, nbit, typeof(V), typeof(VS1),typeof(VT)}(k,t,l,C_tilda,D_tilda,N,V,VS1,VS2,WS1,WS2,VT,WT,F)
    end
end

get_param_type(blk::QuLDEUnitParam{CPType}) where {CPType} = CPType
get_param_type(blk::QuLDEnonUnitParam{CPType}) where {CPType} = CPType

function quldecircuit(blk::QuLDEUnitParam, M::Matrix, nbit::Int, n::Int)
    T = log2i(blk.k+1)
    CPType = get_param_type(blk)
    circuitInit = chain(n, control((-1,),(2:T+1...,)=>blk.VS1),control((1,),(2:T+1...,)=>blk.VS2))

    circuitIntermediate = chain(n)
    a = Array{Int32,1}(undef, T)
    U = Matrix{CPType}(I, 1<<nbit,1<<nbit)

    for i in 0:blk.k
        digits!(a,i,base = 2)
        G = matblock(U)
        push!(circuitIntermediate,control(n, (-1*collect(2:T+1).*((-1*ones(Int, T)).^a)...,), (T+2:n...,)=>G))
        U = M*U
    end

    circuitFinal = chain(n, control((-1,),(2:T+1...,)=>blk.WS1),control((1,),(2:T+1...,)=>blk.WS2), put(1=>blk.V))

    return chain(circuitInit,circuitIntermediate,circuitFinal)
end

function quldecircuit(blk::QuLDEnonUnitParam, M::Matrix, nbit::Int, n::Int)
    circuitL = chain(n)
    for i in 1:blk.k
        push!(circuitL, put(n,(blk.k+2+(i-1)*blk.l:blk.k+2+i*blk.l-1) => blk.VT))
    end

    circuitLinv = chain(n)
    for i in 1:blk.k
        push!(circuitLinv, put(n,(blk.k+2+(i-1)*blk.l:blk.k+2+i*blk.l-1) => blk.WT))
    end

    circuitInit = chain(n, control((-1,), (2:blk.k + 1...,) => blk.VS1), control((1,),(2:blk.k + 1...,)=>blk.VS2), circuitL)

    circuitIntermediate = chain(n)

    a = Array{Int64,1}(undef, blk.l)
    for i in 1:blk.k
        for j in 0:2^blk.l-1
            digits!(a,j,base = 2)
            push!(circuitIntermediate, control(n, (i+1, -1*collect(blk.k+2+(i-1)*blk.l:blk.k+2+i*blk.l-1).*((-1*ones(Int, blk.l)).^a)...,), (n - nbit + 1 : n...,)=>blk.F[j+1]))
        end
    end

    circuitFinal = chain(n, control((-1,), (2:blk.k + 1...,)=>blk.WS1), control((1,), (2:blk.k + 1...,)=>blk.WS2), circuitLinv, put(1=>blk.V))

    return chain(circuitInit, circuitIntermediate, circuitFinal)
end

function DiffEqBase.solve(prob::QuLDEProblem{uType,tType,isinplace, F, P}, alg::QuLDE, k::Int = 3, l::Int = 2) where {uType,tType,isinplace, F, P}
    M = prob.A
    b = prob.b
    t = prob.tspan[2] - prob.tspan[1]
    x = prob.u0
    Mu = M/opnorm(M)
    siz, = size(M)
    nbit = log2i(siz)

    if (isunitary(Mu))
        blk = QuLDEUnitParam(k,t,prob)
        T = log2i(k+1)
        n = 1 + T + nbit
        CPType = get_param_type(blk)
        inreg = ArrayReg(x/norm(x)) ⊗ zero_state(CPType,T) ⊗ ( (blk.C_tilda/blk.N) * zero_state(CPType,1) )+  ArrayReg(b/norm(b)) ⊗ zero_state(CPType,T) ⊗ ((blk.D_tilda/blk.N) * ArrayReg(CPType, bit"1") )
    else
        blk = QuLDEnonUnitParam(k,t,l,prob)
        n = 1 + k*(1 + l) + nbit
        CPType = get_param_type(blk)
        inreg = ArrayReg(x/norm(x))⊗ zero_state(CPType, k*l)  ⊗ zero_state(CPType ,k) ⊗ ( (blk.C_tilda/blk.N) * zero_state(CPType, 1) )+  ArrayReg(b/norm(b)) ⊗ zero_state(CPType, k*l) ⊗ zero_state(CPType, k) ⊗ ((blk.D_tilda/blk.N) * ArrayReg(CPType,bit"1") )
    end
    cir = quldecircuit(blk,Mu,nbit,n)
    res = apply!(inreg,cir) |> focus!(1:n - nbit...,) |> select!(0) |> state

    out = (blk.N^2)*(vec(res))
    return out
end;
