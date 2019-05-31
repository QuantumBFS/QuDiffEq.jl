export array_qudiff, prepare_init_state, bval, aval
import QuAlgorithmZoo.hhlsolve
import LinearAlgebra.eigvals
"""
    Based on : arxiv.org/abs/1010.2745v2

    * array_qudiff(N_t,N,h,A) - generates matrix for k-step solver
    * prepare_init_state(b,x,h,N_t) - generates inital states
    * solve_qudiff - solver

    x' = Ax + b

    * A - input matrix.
    * b - input vector.
    * x - inital vector
    * N - dimension of b (as a power of 2).
    * h - step size.
    * tspan - time span.
"""

"""
    LDEMSAlgHHL
    * step - step for multistep method
    * α - coefficients for xₙ
    * β - coefficent for xₙ'
"""

function bval(g::Function,alg::LDEMSAlgHHL,t,h)
    b = zero(g(1))
    for i in 1:(alg.step)
        b += alg.β[i]*g(t-(i-1)*h)
    end
    return b
end

function aval(g::Function,alg::LDEMSAlgHHL,t,h)
    sz, = size(g(1))
    A = Array{ComplexF64}(undef,sz,(alg.step + 1)*sz)
    i_mat = Matrix{Float64}(I, size(g(1)))
    A[1:sz,sz*(alg.step) + 1:sz*(alg.step + 1)] = i_mat
    for i in 1:alg.step
        A[1:sz,sz*(i - 1) + 1: sz*i] = -1*(alg.α[alg.step - i + 1]*i_mat + h*alg.β[alg.step - i + 1]*g(t - (alg.step - i)*h))
    end
    return A
end

function prepare_init_state(g::Function,alg::LDEMSAlgHHL,tspan::NTuple{2, Float64},x::Vector,h::Float64)
    N_t = round(Int, (tspan[2] - tspan[1])/h + 1) #number of time steps
    N = nextpow(2,2*N_t + 1) # To ensure we have a power of 2 dimension for matrix
    sz, = size(g(1))
    init_state = zeros(ComplexF64,2*(N)*sz)
    #inital value
    init_state[1:sz] = x
    for i in 2:N_t
        b = bval(alg,h*(i - 1) + tspan[1],h) do t g(t) end
        init_state[Int(sz*(i - 1) + 1):Int(sz*(i))] = h*b
    end
    return init_state
end

function array_qudiff(g::Function,alg::LDEMSAlgHHL,tspan::NTuple{2, Float64},h::Float64)
    sz, = size(g(1))
    i_mat = Matrix{Float64}(I, size(g(1)))
    N_t = round(Int, (tspan[2] - tspan[1])/h + 1) #number of time steps
    N = nextpow(2,2*N_t + 1) # To ensure we have a power of 2 dimension for matrix
    A_ = zeros(ComplexF64, N*sz, N*sz)
    # Generates First two rows
    @inbounds A_[1:sz, 1:sz] = i_mat
    @inbounds A_[sz + 1:2*sz, 1:sz] = -1*(i_mat + h*g(tspan[1]))
    @inbounds A_[sz + 1:2*sz,sz+1:sz*2] = i_mat
    #Generates additional rows based on k - step
    for i in 3:alg.step
        @inbounds A_[sz*(i - 1) + 1:sz*i, sz*(i - 3) + 1:sz*i] = aval(QuAB2(),(i-2)*h + tspan[1],h) do t g(t) end
    end
    for i in alg.step + 1:N_t
        @inbounds A_[sz*(i - 1) + 1:sz*(i), sz*(i - alg.step - 1) + 1:sz*i] = aval(alg,(i - 2)*h + tspan[1],h) do t g(t) end
    end
    #Generates half mirroring matrix
    for i in N_t + 1:N
        @inbounds A_[sz*(i - 1) + 1:sz*(i), sz*(i - 2) + 1:sz*(i - 1)] = -1*i_mat
        @inbounds A_[sz*(i - 1) + 1:sz*(i), sz*(i - 1) + 1:sz*i] = i_mat
    end
    A_ = [zero(A_) A_;A_' zero(A_)]
    return A_
end

function DiffEqBase.solve(prob::QuLDEProblem{F,C,U,T}, alg::LDEMSAlgHHL, dt = (prob.tspan[2]-prob.tspan[1])/100, n_reg::Int = 12) where {F,C,U,T}
    A = prob.A
    b = prob.b
    tspan = prob.tspan
    x = prob.u0

    matx = array_qudiff(alg, tspan, dt) do t A(t) end
    initstate = prepare_init_state(alg, tspan, x, dt) do t b(t) end
    λ = maximum(eigvals(matx))
    C_value = minimum(eigvals(matx) .|> abs)*0.01;
    matx = 1/(λ*2)*matx
    initstate = initstate*1/(2*λ) |> normalize!
    res = hhlsolve(matx,initstate, n_reg, C_value)
    res = res/λ
    return res
end;
