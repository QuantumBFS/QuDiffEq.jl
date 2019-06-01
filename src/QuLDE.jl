
struct QuLDEParam{W,L,Q, VType <: GeneralMatrixBlock{1}, SType <: GeneralMatrixBlock{Q}}
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

    function QuLDEParam(k::Int,t::L) where {L}
        C(m) = (t)^(m)/factorial(m)

        C_tilda = 0
        for i in 1:k
            C_tilda = C_tilda + C(i)
        end
        D_tilda = C_tilda
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

function quldecircuit(blk::QuLDEParam, M::Matrix)
    T = Int(log2(blk.k+1))
    siz, = size(M)
    nbit = Int(log2(siz))

    n = 1 + T + nbit
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

function DiffEqBase.solve(prob::QuLDEProblem{F,Q,R,P}, alg::QuLDE, k::Int = 3) where {F,Q,R,P}
    M = prob.A
    b = prob.b
    t = prob.tspan[2]
    x = prob.u0
    T = Int(log2(k+1))

    blk = QuLDEParam(k,t)
    cir = quldecircuit(blk,M)

    inreg = ArrayReg(x) ⊗ zero_state(T) ⊗ ( (blk.C_tilda/blk.N) * zero_state(1) )+  ArrayReg(b) ⊗ zero_state(T) ⊗ ((blk.D_tilda/blk.N) *ArrayReg(bit"1") )

    res = apply!(inreg,cir) |> focus!(1:T+1...,) |> select!(0) |> state

    out = (blk.N^2)*(vec(res))
    return out
end;
