using Yao
using QuDiffEq
using LinearAlgebra
using BitBasis
using OrdinaryDiffEq
using DiffEqOperators
using Plots

# vertex and egde values
vertx = 7
ege = 8

#Constructing incidence matrix B
B = zeros(vertx,ege)
@inbounds for i in 1:vertx
    B[i,i] = -1
    B[i,i+1] = 1
end
B[1,1] = 1

de = 0.5
sn = sin.((0.0:de:(vertx-1)*de)*2*pi/((vertx-1)*de))
u0 = ComplexF64.(sn)

# Intial conditions (stationary to begin with)
u1 = [u0; zero(u0); 0.0; 0.0]
u_ = Float64[u0 zero(u0)]

k = 2   # order in Taylor expansion
t = 1e-2 # time step
B_t = transpose(B)
n = 11  # number of steps
a = 1e-1 # spactial discretization

function make_hamiltonian(B,B_t,a)
    vertx,ege = size(B)
    n = nextpow(2,vertx+ege)
    H = zeros(n,n)
    @inbounds H[1:vertx,vertx+1:vertx+ege] = B
    @inbounds H[vertx+1:vertx+ege,1:vertx] = B_t
    H = -im/a*H
    return H
end

function do_pde(ϕ,B,B_t,k,t,a,n)
    vertx, = size(B)
    H = make_hamiltonian(B,B_t,a)
    res = Array{Array{ComplexF64,1},1}(undef,n)
    res[1] = @view ϕ[1:vertx]
    for i in 2:n
    r, N = taylorsolve(H,ϕ,k,t)
    ϕ = N*vec(state(r))
    res[i] = @view ϕ[1:vertx]
    end
    res_real = real(res)
    return res_real
end

#Dirichlet
res1 = do_pde(u1,B,B_t,k,t,a,n)

#Neumann
#res2 = do_pde(u2,B_t,B,k,t,a,n)

#plot(res1,legend = false)

const Dd = DerivativeOperator{Float64}(2,2,0.1,vertx,:Dirichlet,:Dirichlet)
function f(du,u,p,t)
    buffer, D = p
    u1 = @view(u[:,1])
    u2 = @view(u[:,2])
    mul!(buffer, D, u2)
    Du = buffer

    du[:,1] = Du
    du[:,2] = u1
end

tspan = (0.0,0.1)
prob1 = ODEProblem(f,u_,tspan,(zero(u_[:,1]), Dd))
sol1 = solve(prob1,Tsit5(),dt=0.01,adaptive = false)
s1 = Array{Array{Float64,1},1}(undef,n)
for i in 1:n
    s1[i] = @view (sol1.u[i][:,1])
end

#plot(s1,legend = false)

plot(res - s,legend = false)
