export TaylorParam
export taylorcircuit, get_param_type, taylorsolve, circuit_init, circuit_final, circuit_intermediate
export v,lc

const C(m, x, opn, t, c) = norm(x)*(opn*t*c)^(m)/factorial(m)
const D(m, x, opn, t, c) = norm(x)*(opn*t*c)^(m-1)*t/factorial(m)
struct TaylorParam{CPType, UType, L, HM, VType, S1Type, S2Type, WType}
    k::Int # Taylor expansion upto k
    t::L
    H::HM
    l::Int
    rs::Int
    C_tilda::CPType
    D_tilda::CPType
    N::CPType
    V::VType
    VS1::S1Type # Q is input size of GMblock
    VS2::S2Type
    VT::WType

    function TaylorParam(k::Int, t::L, H::Matrix, x::Array{CPType,1}) where {L,CPType}
        opn = opnorm(H)
        u = isunitary(H/opn)
        VT = nothing
        C_tilda = 0
        if u
            for i in 0:k
                C_tilda = C_tilda + C(i, x, opn,t,1)
            end
            N = C_tilda
            C_tilda = sqrt(C_tilda)
            rs = log2i(k+1)
            VS1 = calc_vs1(k,rs,x,t,opn, C_tilda)
            l = 0
        else
            for i in 0:k
                C_tilda = C_tilda + C(i, x, opn,t,2)
            end
            N = C_tilda
            C_tilda = sqrt(C_tilda)
            rs = k
            l = 2
            VT = calc_vt(CPType)
            VS1 = calc_vs1(k,x,t,opn, C_tilda)
        end

        new{CPType, u, L, Array{CPType,2}, Nothing, typeof(VS1), Nothing, typeof(VT)}(k, t, H, l, rs, C_tilda, zero(C_tilda), N, nothing, VS1, nothing,VT)
    end

    function TaylorParam(k::Int,t::L,prob::QuLDEProblem{uType, tType, isinplace, F, P, T}) where {L,uType, tType, isinplace, F, P, T}
        CPType = eltype(prob.u0)
        opn = opnorm(prob.A)
        u = isunitary(prob.A/opn)
        VS2 = nothing
        V = nothing
        VT = nothing
        D_tilda = 0
        C_tilda = 0

        if u
            for i in 1:k
                D_tilda = D_tilda + D(i, prob.b, opn,t,1)
            end
            N = D_tilda
            D_tilda = sqrt(D_tilda)
            rs = log2i(k+1)
            VS1 = calc_vs2(k,rs,prob.b,t, opn, D_tilda)
        else
            for i in 1:k
                D_tilda = D_tilda + D(i, prob.b, opn,t,2)
            end
            N = D_tilda
            D_tilda = sqrt(D_tilda)
            VT = calc_vt(CPType)
            u = isunitary(prob.A/opn)
            rs = k
            VS1 = calc_vs2(k,prob.b,t, opn, D_tilda)
        end
        if !(T)
            VS2 = VS1
            if u
                for i in 0:k
                    C_tilda = C_tilda + C(i, prob.u0, opn, t,1)
                end
                C_tilda = sqrt(C_tilda)
                VS1 = calc_vs1(k,rs,prob.u0,t, opn, C_tilda)
            else
                for i in 0:k
                    C_tilda = C_tilda + C(i, prob.u0, opn, t,2)
                end
                C_tilda = sqrt(C_tilda)
                VS1 = calc_vs1(k,prob.u0,t, opn, C_tilda)
            end
            N = sqrt(C_tilda^2 + D_tilda^2)
            V = CPType[C_tilda/N D_tilda/N; D_tilda/N -1*C_tilda/N]
            V = matblock(V)
        end
        new{CPType, u, L, Array{CPType,2}, typeof(V), typeof(VS1), typeof(VS2),typeof(VT)}(k, t, prob.A, 2, rs, C_tilda, D_tilda, N, V, VS1, VS2,VT)
    end
end

function calc_vs1(k::Int, x::Vector, t::Real, opn::Real, C_tilda::T) where T
    CPType = eltype(x)
    VS1 = rand(CPType,2^k,2^k)
    @inbounds VS1[:,1] = zero(VS1[:,1])
    @inbounds for j in 0:k
        VS1[(2^k - 2^(k-j) + 1),1] = sqrt(C(j, x, opn, t,2))/C_tilda
    end
    VS1 = -1*qr(VS1).Q
    return VS1
end

function calc_vs1(k::Int, rs::Int, x::Vector, t::Real, opn::Real, C_tilda::T) where T
    CPType = eltype(x)
    @inbounds VS1 = rand(CPType,2^rs,2^rs)
    @inbounds for j in 0:k
        VS1[j+1,1] = sqrt(C(j, x, opn, t,1))/C_tilda
    end
    VS1 = -1*qr(VS1).Q
    return VS1
end

function calc_vs2(k::Int, x::Vector, t::Real, opn::Real, D_tilda::T) where T
    CPType = eltype(x)
    VS2 = rand(CPType,2^k,2^k)
    @inbounds VS2[:,1] = zero(VS2[:,1])
    @inbounds for j in 0:k-1
        VS2[(2^k - 2^(k-j) + 1),1] = sqrt(D(j+1, x, opn, t, 2))/D_tilda
    end
    VS2 = -1*qr(VS2).Q
    return VS2
end

function calc_vs2(k::Int, rs::Int, x::Vector, t::Real, opn::Real, D_tilda::T) where T
    CPType = eltype(x)
    VS2 = rand(CPType,2^rs,2^rs)
    @inbounds for j in 0:k - 1
        VS2[j+1,1] = sqrt(D(j+1, x, opn, t, 1))/D_tilda
    end
    VS2[k+1,1] = 0;
    VS2 = -1*qr(VS2).Q
    return VS2
end

function unitary_decompose(H::Array{T,2}) where T
    if isunitary(H)
        return H
    end
    Mu = H/opnorm(H)
    nbit, = size(H)
    nbit = log2i(nbit)
    B1 = convert(Array{T,2},1/2*(Mu + Mu'))
    B2 = convert(Array{T,2},-im/2*(Mu - Mu'))
    iden = Matrix{T}(I,size(B1))
    F = Array{Array{T,2},1}(undef,4)
    F[1] = B1 + im*sqrt(iden - B1*B1)
    F[2] = B1 - im*sqrt(iden - B1*B1)
    F[3] = im*B2 - sqrt(iden - B2*B2)
    F[4] = im*B2 + sqrt(iden - B2*B2)
    return F

end

function calc_vt(::Type{CPType} = ComplexF32) where CPType
    VT = rand(CPType,4,4)
    VT[:,1] = 0.5*ones(4)
    VT = -1*qr(VT).Q
    return VT
end

v(n::Int,c::Int, T::Int, V::AbstractMatrix) = put(n, (c + 1:T + c...,)=>matblock(V))
v(n::Int,c::Int, j::Tuple, T::Int, V::AbstractMatrix) = control(n, j,(1 + c:T + c...,)=>matblock(V))
lc(n::Int,c::Int, i::Int, k::Int,l::Int, V::AbstractMatrix) = put(n, (k+c+1+(i-1)*l:k+1+c+i*l-1) => matblock(V))

circuit_init(n::Int, c::Int, blk::TaylorParam{CPType, true}) where CPType = chain(n, v(n, c,(-1,), blk.rs, blk.VS1),v(n, c, (1,), blk.rs, blk.VS2))
circuit_init(n::Int, blk::TaylorParam{CPType, true}) where CPType = chain(n, v(n, 0, blk.rs, blk.VS1))

function circuit_init(n::Int,c::Int, blk::TaylorParam{CPType, false}) where CPType
    cir = chain(n, v(n, c,(-1,), blk.rs, blk.VS1),v(n, c, (1,), blk.rs, blk.VS2))
    for i in 1:blk.k
        push!(cir, lc(n, c, i, blk.k, blk.l, blk.VT))
    end
    return cir
end

function circuit_init(n::Int, blk::TaylorParam{CPType, false}) where CPType
    cir = chain(n, v(n, 0, blk.rs, blk.VS1))
    for i in 1:blk.k
        push!(cir, lc(n, 0, i, blk.k, blk.l, blk.VT))
    end
    return cir
end

function circuit_intermediate(n::Int, c::Int, blk::TaylorParam{CPType, true}) where CPType
    H = blk.H/opnorm(blk.H)
    k = blk.k
    l = blk.l
    rs = blk.rs
    cir = chain(n)
    nbit = n - rs - c
    a = Array{Int32,1}(undef, rs)
    U = Matrix{CPType}(I, 1<<nbit,1<<nbit)

    for i in 0:k
        digits!(a,i,base = 2)
        G = matblock(U)
        push!(cir,control(n, (-1*collect(1+c:rs+c).*((-1*ones(Int, rs)).^a)...,), (rs+1+c:n...,)=>G))
        U = H*U
    end
    return cir
end

function circuit_intermediate(n::Int, c::Int, blk::TaylorParam{CPType, false}) where CPType
    H = blk.H/opnorm(blk.H)
    k = blk.k
    l = blk.l
    rs = blk.rs
    cir = chain(n)
    nbit = n - k*(l+1) - c
    F = unitary_decompose(H)
    a = Array{Int64,1}(undef, l)
    for i in 1:k
        for j in 0:2^l-1
            digits!(a,j,base = 2)
            push!(cir, control(n, (i + c, -1*collect(k+1+c+(i-1)*l:k+1+c+i*l-1).*((-1*ones(Int, l)).^a)...,), (n - nbit + 1 : n...,)=>matblock(F[j+1])))
        end
    end
    return cir
end

circuit_final(n::Int, c::Int, blk::TaylorParam{CPType, true}) where CPType = chain(n, v(n, c, (-1,), blk.rs, blk.VS1'),v(n, c, (1,),blk.rs, blk.VS2'))
circuit_final(n::Int, blk::TaylorParam{CPType, true}) where CPType = chain(n, v(n, 0, blk.rs,blk.VS1'))

function circuit_final(n::Int, c::Int, blk::TaylorParam{CPType, false}) where CPType
    cir = chain(n, v(n, c, (-1,), blk.rs, blk.VS1'),v(n, c, (1,),blk.rs, blk.VS2'))
    for i in 1:blk.k
        push!(cir, lc(n,c,i,blk.k,blk.l,blk.VT'))
    end
    return cir
end

function circuit_final(n::Int, blk::TaylorParam{CPType, false}) where CPType
    cir = chain(n, v(n, 0, blk.rs,blk.VS1'))
    for i in 1:blk.k
        push!(cir, lc(n,0,i,blk.k,blk.l,blk.VT'))
    end
    return cir
end

get_param_type(blk::TaylorParam{CPType}) where {CPType} = CPType

function taylorcircuit(n::Int, blk::TaylorParam)
    circinit = circuit_init(n,blk)
    circmid = circuit_intermediate(n,0,blk)
    circfin =circuit_final(n,blk)
    return chain(circinit, circmid, circfin)
end

function taylorsolve(H::Matrix, x::Vector, k::Int, t::Real)

    nbit = log2i(length(x))
    blk = TaylorParam(k,t,H,x)
    CPType = get_param_type(blk)
    rs = blk.rs
    k = blk.k
    l = blk.l
    if rs != k
        n = rs + nbit
        inreg = ArrayReg(x/norm(x)) ⊗ zero_state(CPType,rs)
    else
        n = k*(1 + l) + nbit
        CPType = get_param_type(blk)
        inreg = ArrayReg(x/norm(x))⊗ zero_state(CPType, k*l)  ⊗ zero_state(CPType, k)
    end
    cir = taylorcircuit(n, blk)
    r = apply!(inreg,cir) |> focus!(1:n - nbit...,) |> select!(0)
    return r, blk.N
end
