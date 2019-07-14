export func_transform, euler_matrix, nonlinear_transform, make_input_vector, make_unitary

function make_input_vector(x::Vector)
    if !(norm(x) ≈ 1)
        throw(ArgumentError("Input vector is not normalized"))
    end
    CPType = eltype(x)
    len = length(x)
    siz = nextpow(2,len + 1)
    z = zeros(CPType, siz)
    z[1] = CPType(1)
    z[2:len+1] = x
    normalize!(z)
    reg = [1,0] ⊗ z ⊗ z
    return reg
end

make_unitary(A::Matrix) = [zero(A) im*A'; -im*A zero(A)]

function nonlinear_transform(H::Matrix, x::Vector, k::Int, ϵ::Real = 1e-4)
    r, N = taylorsolve(im*H,x,k,ϵ)
    n = log2i(length(x))
    nb = Int((n-1)/2)
    r = relax!(r) |> focus!(n)|> select!(1) |> focus!(1:nb...,) |> select!(0)
    return r,sqrt(2)*N/ϵ
end

function func_transform(A::Matrix, x::Vector, k::Int,ϵ::Real = 1e-4)
    reg = make_input_vector(x)
    H = make_unitary(A)
    r, N = nonlinear_transform(H,reg,k,ϵ)
    return r, N
end
function euler_matrix(A::Matrix,b::Vector,h::Real)
    n = length(b)
    A = h*A
    A[1,1] = 1
    for i in 1:n
        A[4*i+ 1, i+1] += 1
    end
    return A
end

function euler_matrix_update(A::Matrix,b::Vector,nrm::Real)
    n = length(b)
    A = nrm^2*A
    A[1,1] = 1
    for i in 1:n
        A[4*i+ 1, 1] = A[4*i+ 1, 1]/(nrm^2)
        @. A[4*i+ 1, 2:n+1] = A[4*i+ 1, 2:n+1]/nrm
        for j in 1:n
            A[4*i + 1, 4*j + 1] = A[4*i + 1, 4*j + 1]/nrm
        end
    end
    return A
end

function DiffEqBase.solve(prob::QuODEProblem,alg::QuNLDE; dt = (prob.tspan[2]-prob.tspan[1])/10, kwargs...)
    A = prob.A
    b = prob.b
    k = alg.k
    ϵ = alg.ϵ
    tspan = prob.tspan
    len = round(Int,(tspan[2] - tspan[1])/dt) + 1
    siz = length(b)
    res = Array{eltype(b),2}(undef,len,siz)
    res[1,:] = b
    A = euler_matrix(A,b,dt)
    H = make_unitary(A)
    reg = make_input_vector(b)
    r, N = nonlinear_transform(H,reg,k,ϵ)
    ntem = 1
    for step in 2:len - 1
        tem = vec(state(r))*N*sqrt(2)
        res[step,:] = tem[2:siz+1]
        ntem = norm(res[step,:])
        C = euler_matrix_update(A,b,ntem)
        H = make_unitary(C)
        reg = res[step,:]/ntem
        reg = make_input_vector(reg)
        r, N = nonlinear_transform(H,reg,k,ϵ)
    end
    tem = vec(state(r))*N*sqrt(2)
    res[len,:] = tem[2:siz+1]
    return res
end
