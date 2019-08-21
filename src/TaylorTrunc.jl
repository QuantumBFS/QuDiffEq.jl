export TaylorParam
export taylorcircuit, taylorsolve, circuit_ends, circuit_intermediate
export v,lc
#export get_param_type, circuit_final

const C(m, x, opn, t, c) = norm(x)*(opn*t*c)^(m)/factorial(m)
const D(m, x, opn, t, c) = norm(x)*(opn*t*c)^(m-1)*t/factorial(m)
"""

    TaylorParam(k::Int, t::L, H::Matrix, x::Array{CPType,1})

Assigns values to parameters required for Taylor series based Hamiltonian simulation.

    * k : sets order of Taylor series expansion
    * t : time of evolution
    * H : Hamiltonian
    * x : intial vector

    TaylorParam(k::Int,t::L,prob::QuLDEProblem{uType, tType, isinplace, F, P, T})

    * k : sets order of Taylor series expansion
    * t : time of evolution
    * prob : wrapper for Linear differential equation problems (contains H, x, b)
"""
struct TaylorParam{CPType, UType, L, HM}
    k::Int # Taylor expansion upto k
    t::L
    H::HM
    l::Int
    rs::Int
    C_tilda::CPType
    D_tilda::CPType
    N::CPType
    function TaylorParam(k::Int, t::L, H::Matrix, x::Array{CPType,1}) where {L,CPType}
        opn = opnorm(H)
        u = isunitary(H/opn)
        C_tilda = 0
        if u
            c = 1
            rs = log2i(k+1)
            l = 0
        else
            c = 2
            rs = k
            l = 2
        end

        for i in 0:k
            C_tilda = C_tilda + C(i, x, opn,t,c)
        end
        N = C_tilda
        C_tilda = sqrt(C_tilda)
        new{CPType, u, L, Array{CPType,2}}(k, t, H, l, rs, C_tilda, zero(C_tilda), N)
    end

    function TaylorParam(k::Int,t::L,prob::QuLDEProblem{uType, tType, isinplace, F, P, T}) where {L,uType, tType, isinplace, F, P, T}
        CPType = eltype(prob.A)
        opn = opnorm(prob.A)
        u = isunitary(prob.A/opn)
        D_tilda = 0
        C_tilda = 0
        if u
            c = 1
            rs = log2i(k+1)
            l = 0
        else
            c = 2
            rs = k
            l = 2
        end
        for i in 1:k
            D_tilda = D_tilda + D(i, prob.b, opn,t,c)
        end
        N = D_tilda
        D_tilda = sqrt(D_tilda)
        if !(T)
            for i in 0:k
                C_tilda = C_tilda + C(i, prob.u0, opn, t,c)
            end
            C_tilda = sqrt(C_tilda)
            N = sqrt(C_tilda^2 + D_tilda^2)
        end
        new{CPType, u, L, Array{CPType,2}}(k, t, prob.A, l, rs, C_tilda, D_tilda, N)
    end
end

"""
    calc_vs1(blk::TaylorParam, x::Vector{CPType}, opn::Real) ->  Matrix{CPType}

Calculates VS1 block for Taylor circuit.
"""
function calc_vs1(blk::TaylorParam, x::Vector{CPType}, opn::Real) where CPType
    k = blk.k
    rs = blk.rs
    t = blk.t
    C_tilda = blk.C_tilda
    VS1 = rand(CPType,2^rs,2^rs)
    @inbounds VS1[:,1] = zero(VS1[:,1])
    if rs == k
        @inbounds for j in 0:k
            VS1[(2^k - 2^(k-j) + 1),1] = sqrt(C(j, x, opn, t,2))/C_tilda
        end
    else
        @inbounds for j in 0:k
            VS1[j+1,1] = sqrt(C(j, x, opn, t,1))/C_tilda
        end
    end
    VS1 = -1*qr(VS1).Q
    return VS1
end

"""
    calc_vs2(blk::TaylorParam, x::Vector{CPType}, opn::Real) ->  Matrix{CPType}

    opn: operator norm of the input matrix.

Calculates VS2 block for Taylor circuit.
"""
function calc_vs2(blk::TaylorParam, x::Vector{CPType}, opn::Real) where CPType
    k = blk.k
    rs = blk.rs
    t = blk.t
    D_tilda = blk.D_tilda
    VS2 = rand(CPType,2^rs,2^rs)
    @inbounds VS2[:,1] = zero(VS2[:,1])
    if rs == k
        @inbounds for j in 0:k-1
            VS2[(2^k - 2^(k-j) + 1),1] = sqrt(D(j+1, x, opn, t, 2))/D_tilda
        end
    else
        @inbounds for j in 0:k - 1
            VS2[j+1,1] = sqrt(D(j+1, x, opn, t, 1))/D_tilda
        end
    end
    VS2 = -1*qr(VS2).Q
    return VS2
end

"""
    unitary_decompose(H::Array{T,2}) -> Array{Array{T,2},1}

Generates a linear compostion of unitary matrices for argument H.
"""
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

"""
    calc_vt(T)

    calc_vt(::Type{CPType}) -> Matrix{CPType}

Generates VT block for non-unitary Taylor circuit.
"""
function calc_vt(::Type{CPType} = ComplexF32) where CPType
    VT = rand(CPType,4,4)
    VT[:,1] = 0.5*ones(4)
    VT = -1*qr(VT).Q
    return VT
end

"""
    v(n::Int,c::Int, T::Int, V::AbstractMatrix) -> ChainBlock{n}

Builds T input block.
n : total number of qubits
V : T input matrix
c : starting qubit
"""
v(n::Int,c::Int, T::Int, V::AbstractMatrix) = concentrate(n, matblock(V), (c + 1:T + c...,))

"""
    v(n::Int,c::Int, j::Tuple, T::Int, V::AbstractMatrix) -> ControlBlock{n}

Builds a T input, j - control block.
n : total number of qubits
V : T input matrix
c : starting qubit
j : control bits tuple
"""
v(n::Int,c::Int, j::Tuple, T::Int, V::AbstractMatrix) = control(n, j,(1 + c:T + c...,)=>matblock(V))
lc(n::Int,c::Int, i::Int, k::Int,l::Int, V::AbstractMatrix) = concentrate(n, matblock(V), (k+c+1+(i-1)*l:k+1+c+i*l-1...,))


"""
    circuit_ends(n::Int, blk::TaylorParam{CPType, true}, VS1::AbstractMatrix, VS2::AbstractMatrix)

Generates the part of circuit that computes and decomputes the superposition of the ancilla bits, in unitary H `quldecircuit`
"""
circuit_ends(n::Int, blk::TaylorParam{CPType, true}, VS1::AbstractMatrix, VS2::AbstractMatrix) where CPType = chain(n, v(n, 1,(-1,), blk.rs, VS1),v(n, 1, (1,), blk.rs, VS2))


"""
    circuit_ends(n::Int, blk::TaylorParam{CPType, true}, VS1::AbstractMatrix)

Generates the part of circuit that computes and decomputes the superposition of the ancilla bits, in unitary H `taylorcircuit`
"""
circuit_ends(n::Int, blk::TaylorParam{CPType, true}, VS1::AbstractMatrix) where CPType = chain(n, v(n, 0, blk.rs, VS1))

"""
    circuit_ends(n::Int, blk::TaylorParam{CPType, false}, VS1::AbstractMatrix, VS2::AbstractMatrix, VT::AbstractMatrix)

Generates the part of circuit that computes and decomputes the superposition of the ancilla bits, in non-unitary H `quldecircuit`
"""
function circuit_ends(n::Int, blk::TaylorParam{CPType, false}, VS1::AbstractMatrix, VS2::AbstractMatrix, VT::AbstractMatrix) where CPType
    cir = chain(n, v(n, 1, (-1,), blk.rs, VS1),v(n, 1, (1,),blk.rs, VS2))
    for i in 1:blk.k
        push!(cir, lc(n,1,i,blk.k,blk.l,VT))
    end
    return cir
end


"""
    circuit_ends(n::Int, blk::TaylorParam{CPType, false}, VS1::AbstractMatrix, VT::AbstractMatrix)

Generates the part of circuit that computes and decomputes the superposition of the ancilla bits, in non-unitary H `taylorcircuit`
"""
function circuit_ends(n::Int, blk::TaylorParam{CPType, false}, VS1::AbstractMatrix, VT::AbstractMatrix) where CPType
    cir = chain(n, v(n, 0, blk.rs,VS1))
    for i in 1:blk.k
        push!(cir, lc(n,0,i,blk.k,blk.l,VT))
    end
    return cir
end
"""
    circuit_intermediate(n::Int, c::Int, blk::TaylorParam{CPType, true})

Generates the intermediate part of the circuit for unitary H.
"""
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
"""
    circuit_intermediate(n::Int, c::Int, blk::TaylorParam{CPType, false})

Generates the intermediate part of the circuit for non-unitary H.
"""
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

"""
    taylorcircuit(n::Int, blk::TaylorParam, VS1::Matrix) ->  ChainBlock{n}

Generates circuit for a unitary H input.
"""
function taylorcircuit(n::Int, blk::TaylorParam, VS1::Matrix)
    circinit = circuit_ends(n,blk,VS1)
    circmid = circuit_intermediate(n,0,blk)
    circfin =circuit_ends(n,blk,VS1')
    return chain(circinit, circmid, circfin)
end
"""
    taylorcircuit(n::Int, blk::TaylorParam, VS1::Matrix, VT::Matrix) ->->  ChainBlock{n}

Generates circuit for a non-unitary H input.
"""
function taylorcircuit(n::Int, blk::TaylorParam, VS1::Matrix, VT::Matrix)
    circinit = circuit_ends(n,blk,VS1,VT)
    circmid = circuit_intermediate(n,0,blk)
    circfin =circuit_ends(n,blk,VS1',VT')
    return chain(circinit, circmid, circfin)
end

"""
    taylorsolve(H::Array{CPType,2}, x::Vector{CPType}, k::Int, t::Real) -> ArrayReg, CPType

Simulates a Hamiltonian using the Taylor truncation method. Returns the state register and inverse probability of finding it.
"""
function taylorsolve(H::Array{CPType,2}, x::Vector{CPType}, k::Int, t::Real) where CPType
    opn = opnorm(H)
    nbit = log2i(length(x))
    blk = TaylorParam(k,t,H,x)
    VS1 = calc_vs1(blk,x,opn)
    rs = blk.rs
    k = blk.k
    l = blk.l
    if rs != k
        n = rs + nbit
        inreg = ArrayReg(x/norm(x)) ⊗ zero_state(CPType,rs)
        cir = taylorcircuit(n, blk, VS1)
    else
        n = k*(1 + l) + nbit
        VT = calc_vt(CPType)
        inreg = ArrayReg(x/norm(x)) ⊗ zero_state(CPType, k*l)  ⊗ zero_state(CPType, k)
        cir = taylorcircuit(n, blk, VS1, VT)
    end
    r = apply!(inreg,cir) |> focus!(1:n - nbit...,) |> select!(0)
    return r, blk.N
end
