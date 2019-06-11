
struct QuLDEUnitParam{W,L,Q, VType <: GeneralMatrixBlock{1}, SType <: GeneralMatrixBlock{Q}}
    k::Int
    t::L
    C_tilda::W
    D_tilda::W
    N::W
    V::VType
    VS1::SType
    VS2::SType
    WS1::SType
    WS2::SType

    function QuLDEUnitParam(k::Int,t::L,prob::QuLDEProblem) where {L}
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
        V = [C_tilda D_tilda; D_tilda -1*C_tilda]/N
        V = convert(Array{ComplexF64,2},V)
        V = matblock(V)

        col1 = (1/C_tilda) * [sqrt(C(m)) for m in 0:k]
        col2 = (1/D_tilda) * [m < k+1 ? sqrt(C(m)) : 0 for m in 1:k+1]
        VS1 = rand(ComplexF64,k+1,k+1)
        VS2 = rand(ComplexF64,k+1,k+1)
        VS1[:,1] = col1
        VS2[:,1] = col2
        VS1 = -1*qr(VS1).Q
        VS2 = -1*qr(VS2).Q
        WS1 = convert(typeof(VS1),VS1')
        WS2 = convert(typeof(VS1),VS2')
        VS1 = matblock(VS1)
        VS2 = matblock(VS2)
        WS1 = matblock(WS1)
        WS2 = matblock(WS2)
        n = Int(log2(k+1))
        new{ComplexF64,L,n, typeof(V), typeof(VS1)}(k,t,C_tilda,D_tilda,N,V,VS1,VS2,WS1,WS2)
    end
end

struct QuLDEnonUnitParam{W,L,Q,NL,Nbit, VType <: GeneralMatrixBlock{1}, SType <: GeneralMatrixBlock{Q}, LType <: GeneralMatrixBlock{NL}}
    k::Int
    t::L
    l::Int
    C_tilda::W
    D_tilda::W
    N::W
    V::VType
    VS1::SType
    VS2::SType
    WS1::SType
    WS2::SType
    VT::LType
    WT::LType
    F::Array{GeneralMatrixBlock{Nbit,Nbit,Complex{Float64},Array{Complex{Float64},2}},1}

    function QuLDEnonUnitParam(k::Int,t::L,l::Int,prob::QuLDEProblem) where {L}
        C(m) = norm(prob.u0)*((opnorm(prob.A)*t*2)^(m))/factorial(m)
        D(n) = norm(prob.b)*(opnorm(prob.A)*t*2)^(n-1)*t/factorial(n)
        nbit, = size(prob.u0)
        nbit = Int(log2(nbit))
        Mu = prob.A/opnorm(prob.A)
        B1 = 1/2*(Mu + Mu')
        B2 = -im/2*(Mu - Mu')
        iden = Matrix{ComplexF64}(I,size(B1))
        F = Array{GeneralMatrixBlock{nbit,nbit,Complex{Float64},Array{Complex{Float64},2}},1}(undef,4)
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
        N = sqrt(C_tilda^2 + D_tilda^2)
        V = [C_tilda D_tilda; D_tilda -1*C_tilda]/N
        V = convert(Array{ComplexF64,2},V)
        V = matblock(V)
        #tested

        col1 = zeros(2^k,1)
        col2 = zeros(2^k,1)
        temp = map(0:k) do j 2^k - 2^(k-j) + 1 end
        cont = 0
        for i in temp
            col1[i] = sqrt(C(cont))
            cont += 1
        end
        cont = 1
        for i in temp[2:end]
            col2[i] = sqrt(D(cont))
            cont += 1
        end
        #tested

        VS1 = rand(ComplexF64,2^k,2^k)
        VS2 = rand(ComplexF64,2^k,2^k)
        VS1[:,1] = col1/sqrt(sum(col1.*col1))
        VS2[:,1] = col2/sqrt(sum(col2.*col2))
        VS1 = -1*qr(VS1).Q
        VS2 = -1*qr(VS2).Q
        WS1 = convert(typeof(VS1),VS1')
        WS2 = convert(typeof(VS1),VS2')
        VS1 = matblock(VS1)
        VS2 = matblock(VS2)
        WS1 = matblock(WS1)
        WS2 = matblock(WS2)

        vcol = [sqrt(0.5) for i in 1:2^l]
        VT = rand(ComplexF64,2^l,2^l)
        VT[:,1] = vcol/sqrt(sum(vcol.*vcol))
        VT = -1*qr(VT).Q
        WT = convert(typeof(VT),VT')
        VT = matblock(VT)
        WT = matblock(WT)
        new{ComplexF64,L,k,l,nbit, typeof(V), typeof(VS1),typeof(VT)}(k,t,l,C_tilda,D_tilda,N,V,VS1,VS2,WS1,WS2,VT,WT,F)
    end
end

function quldecircuit(blk::QuLDEUnitParam, M::Matrix, nbit::Int, n::Int)
    T = Int(log2(blk.k+1))
    circuitInit = chain(n, control((-1,),(2:T+1...,)=>blk.VS1),control((1,),(2:T+1...,)=>blk.VS2))

    circuitIntermediate = chain(n)
    a = Array{Int64,1}(undef, T)
    U = Matrix{ComplexF64}(I, 1<<nbit,1<<nbit)

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
    U = Matrix{ComplexF64}(I, 1<<nbit, 1<<nbit)

    for i in 1:blk.k
        for j in 0:2^blk.l-1
            digits!(a,j,base = 2)
            push!(circuitIntermediate, control(n, (i+1, -1*collect(blk.k+2+(i-1)*blk.l:blk.k+2+i*blk.l-1).*((-1*ones(Int, blk.l)).^a)...,), (n - nbit + 1 : n...,)=>blk.F[j+1]))
        end
    end

    circuitFinal = chain(n, control((-1,), (2:blk.k + 1...,)=>blk.WS1), control((1,), (2:blk.k + 1...,)=>blk.WS2), circuitLinv, put(1=>blk.V))
    return chain(circuitInit, circuitIntermediate, circuitFinal)
end

function DiffEqBase.solve(prob::QuLDEProblem{F,Q,R,P}, alg::QuLDE, k::Int = 3, l::Int = 2) where {F,Q,R,P}
    M = prob.A
    b = prob.b
    t = prob.tspan[2]
    x = prob.u0
    Mu = M/opnorm(M)
    siz, = size(M)
    nbit = Int(log2(siz))

    if (isunitary(Mu))
        blk = QuLDEUnitParam(k,t,prob)
        T = Int(log2(k+1))
        n = 1 + T + nbit
        inreg = ArrayReg(x/norm(x)) ⊗ zero_state(T) ⊗ ( (blk.C_tilda/blk.N) * zero_state(1) )+  ArrayReg(b/norm(b)) ⊗ zero_state(T) ⊗ ((blk.D_tilda/blk.N) *ArrayReg(bit"1") )
    else
        blk = QuLDEnonUnitParam(k,t,l,prob)
        n = 1 + k*(1 + l) + nbit
        inreg = ArrayReg(x/norm(x))⊗ zero_state(k*l)  ⊗ zero_state(k) ⊗ ( (blk.C_tilda/blk.N) * zero_state(1) )+  ArrayReg(b/norm(b)) ⊗ zero_state(k*l) ⊗ zero_state(k) ⊗ ((blk.D_tilda/blk.N) * ArrayReg(bit"1") )
    end
    cir = quldecircuit(blk,Mu,nbit,n)
    res = apply!(inreg,cir) |> focus!(1:n - nbit...,) |> select!(0) |> state

    out = (blk.N^2)*(vec(res))
    return out
end;
