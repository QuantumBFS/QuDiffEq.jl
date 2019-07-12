using Yao
using QuDiffEq
using LinearAlgebra
using BitBasis
using Random
using Test
using YaoBlocks
using OrdinaryDiffEq

#Linear Diff Equation Unitary M
function f(du,u,p,t)
    du[1] = -3*u[1]^2 + u[2]
    du[2] = -u[2]^2 - u[1]*u[2]
end
@testset "QuNLDE_Test" begin
    Random.seed!(4)
    N = 2
    k = 3
    siz = nextpow(2, N + 1)
    x = normalize!(rand(N))

    A = zeros(ComplexF32,2^(siz),2^(siz))
    A[1,1] = ComplexF32(1)
    A[5,3] = ComplexF32(1)
    A[5,6] = ComplexF32(-3)
    A[9,11] = ComplexF32(-1)
    A[9,7] = ComplexF32(-1)
    tspan = (0.0,0.4)
    qprob = QuLDEProblem(A, x, tspan)
    r,N = func_transform(qprob.A, qprob.b, k)
    out = N*vec(state(r))
    r_out = zero(x)
    f(r_out, x,1,1)
    @test isapprox.(r_out, out[2:3]*sqrt(2), atol = 1e-3) |> all

    prob = ODEProblem(f, x, tspan)
    sol = solve(prob, Euler(), dt = 0.1, adaptive = false)
    r_out = transpose(hcat(sol.u...))
    out = solve(qprob, QuNLDE(5), dt = 0.1)
    @test_broken isapprox.(r_out,real(out), atol = 1e-3) |> all
end
