export TaylorUnitParam, TaylornonUnitParam
export taylorcircuit, get_param_type, taylorsolve
struct TaylorUnitParam{CPType, L, HM, Q, SType <: GeneralMatrixBlock{Q}}
    k::Int # Taylor expansion upto k
    t::L
    H::HM # time
    C_tilda::CPType
    VS1::SType # Q is input size of GMblock
    WS1::SType # Q

    function TaylorUnitParam(k::Int, t::L, H::Matrix, x::Array{CPType,1}) where {L,CPType}
        C(m) = norm(x)*(opnorm(H)*t)^(m)/factorial(m)
        C_tilda = 0
        for i in 0:k
            C_tilda = C_tilda + C(i)
        end

        VS1 = rand(CPType,k+1,k+1)
        for j in 0:k
            VS1[j+1,1] = sqrt(C(j)/C_tilda)
        end
        VS1 = -1*qr(VS1).Q
        WS1 = convert(typeof(VS1),VS1')
        VS1 = matblock(VS1)
        WS1 = matblock(WS1)

        n = log2i(k+1)
        new{CPType, L, Array{CPType,2}, n, typeof(VS1)}(k, t, H, C_tilda, VS1, WS1)
    end
end
struct TaylornonUnitParam{CPType, L, HM, Q, NL, Nbit, SType <: GeneralMatrixBlock{Q}, LType <: GeneralMatrixBlock{NL}}
    k::Int # Taylor expansion upto k
    t::L # time
    l::Int
    H::HM # input size of LType block set to 2
    C_tilda::CPType #G1 as per the article
    VS1::SType # Q is input size of GMblock
    WS1::SType # Q
    VT::LType # NL is input size of GMblock
    WT::LType # NL
    F::Array{GeneralMatrixBlock{Nbit,Nbit,CPType,Array{CPType,2}},1} # Nbit is input size

    function TaylornonUnitParam(k::Int,t::L,l::Int,H::Matrix,x::Array{CPType,1}) where {L,CPType}
        C(m) = norm(x)*((opnorm(H)*t*2)^(m))/factorial(m) # alphas equal 1/2
        nbit, = size(x)
        nbit = log2i(nbit)
        Mu = H/opnorm(H)
        B1 = convert(Array{CPType,2},1/2*(Mu + Mu'))
        B2 = convert(Array{CPType,2},-im/2*(Mu - Mu'))
        iden = Matrix{CPType}(I,size(B1))
        F = Array{GeneralMatrixBlock{nbit,nbit,CPType,Array{CPType,2}},1}(undef,4)
        F[1] = matblock(B1 + im*sqrt(iden - B1*B1))
        F[2] = matblock(B1 - im*sqrt(iden - B1*B1))
        F[3] = matblock(im*B2 - sqrt(iden - B2*B2))
        F[4] = matblock(im*B2 + sqrt(iden - B2*B2))

        C_tilda = 0
        for i in 0:k
            C_tilda = C_tilda + C(i)
        end

        VS1 = rand(CPType,2^k,2^k)
        VS1[:,1] = zero(VS1[:,1])
        for j in 0:k
            VS1[(2^k - 2^(k-j) + 1),1] = sqrt(C(j)/C_tilda)
        end
        VS1 = -1*qr(VS1).Q
        WS1 = convert(typeof(VS1),VS1')
        VS1 = matblock(VS1)
        WS1 = matblock(WS1)

        VT = rand(CPType,2^l,2^l)
        VT[:,1] = 0.5*ones(2^l,1)
        VT = -1*qr(VT).Q
        WT = convert(typeof(VT),VT')
        VT = matblock(VT)
        WT = matblock(WT)
        new{CPType, L, Array{CPType,2}, k, l, nbit, typeof(VS1),typeof(VT)}(k, t, l, H, C_tilda, VS1, WS1, VT, WT, F)
    end
end

get_param_type(blk::TaylorUnitParam{CPType}) where {CPType} = CPType
get_param_type(blk::TaylornonUnitParam{CPType}) where {CPType} = CPType

function taylorcircuit(blk::TaylorUnitParam, nbit::Int, n::Int)
    T = log2i(blk.k+1)
    CPType = get_param_type(blk)
    circuitInit = chain(n, put((1:T...,)=>blk.VS1))

    circuitIntermediate = chain(n)
    a = Array{Int32,1}(undef, T)
    U = Matrix{CPType}(I, 1<<nbit,1<<nbit)
    H = blk.H/opnorm(blk.H)
    for i in 0:blk.k
        digits!(a,i,base = 2)
        G = matblock(U)
        push!(circuitIntermediate,control(n, (-1*collect(1:T).*((-1*ones(Int, T)).^a)...,), (T+1:n...,)=>G))
        U = H*U
    end

    circuitFinal = chain(n, put((1:T...,)=>blk.WS1))

    return chain(circuitInit,circuitIntermediate,circuitFinal)
end

function taylorcircuit(blk::TaylornonUnitParam, nbit::Int, n::Int)
    circuitL = chain(n)
    for i in 1:blk.k
        push!(circuitL, put(n,(blk.k+1+(i-1)*blk.l:blk.k+1+i*blk.l-1) => blk.VT))
    end

    circuitLinv = chain(n)
    for i in 1:blk.k
        push!(circuitLinv, put(n,(blk.k+1+(i-1)*blk.l:blk.k+1+i*blk.l-1) => blk.WT))
    end

    circuitInit = chain(n, put((1:blk.k...,) => blk.VS1), circuitL)

    circuitIntermediate = chain(n)

    a = Array{Int64,1}(undef, blk.l)
    for i in 1:blk.k
        for j in 0:2^blk.l-1
            digits!(a,j,base = 2)
            push!(circuitIntermediate, control(n, (i, -1*collect(blk.k+1+(i-1)*blk.l:blk.k+1+i*blk.l-1).*((-1*ones(Int, blk.l)).^a)...,), (n - nbit + 1 : n...,)=>blk.F[j+1]))
        end
    end

    circuitFinal = chain(n, put((1:blk.k...,)=>blk.WS1), circuitLinv)

    return chain(circuitInit, circuitIntermediate, circuitFinal)
end

function taylorsolve(H::Matrix, x::Vector, k::Int, t::Real)
    Mu = H/opnorm(H)
    nbit = log2i(length(x))
    if (isunitary(H))
        blk = TaylorUnitParam(k,t,H,x)
        T = log2i(k+1)
        n = T + nbit
        CPType = get_param_type(blk)
        inreg = ArrayReg(x/norm(x)) ⊗ zero_state(CPType,T)
    else
        l = 2
        blk = TaylornonUnitParam(k,t,l,H,x)
        n = k*(1 + l) + nbit
        CPType = get_param_type(blk)
        inreg = ArrayReg(x/norm(x))⊗ zero_state(CPType, k*l)  ⊗ zero_state(CPType, k)
    end
    cir = taylorcircuit(blk,nbit,n)
    r = apply!(inreg,cir) |> focus!(1:n - nbit...,) |> select!(0)
    return r, blk.C_tilda
end
