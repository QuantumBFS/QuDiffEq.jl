using Yao
using Yao.BitBasis
using Random
using Test, LinearAlgebra
using OrdinaryDiffEq
using QuDiffEq

function diffeq_problem(nbit::Int)
    siz = 1<<nbit
    A = (rand(ComplexF64, siz,siz))
    A = (A + A')/2
    b = normalize!(rand(ComplexF64, siz))
    x = normalize!(rand(ComplexF64, siz))
    A, b, x
end

@testset "Linear_differential_equations_HHL" begin
    Random.seed!(2)
    N = 1
    h = 0.1
    tspan = (0.0,0.6)
    N_t = round(Int, 2*(tspan[2] - tspan[1])/h + 3)
    M, v, x = diffeq_problem(N)
    A(t) = M
    b(t) = v
    nreg = 12
    f(u,p,t) = M*u + v;
    prob = ODEProblem(f, x, tspan)
    qprob = QuLDEProblem(A, b, x, tspan)
    sol = solve(prob, Tsit5(), dt = h, adaptive = false)
    s = vcat(sol.u...)

    res = solve(qprob, QuEuler(nreg), dt = h)
    @test isapprox.(s, res, atol = 0.5) |> all

    res = solve(qprob, QuLeapfrog(nreg), dt = h)
    @test isapprox.(s, res, atol = 0.3) |> all

    res = solve(qprob, QuAB2(nreg), dt = h)
    @test isapprox.(s, res, atol = 0.3) |> all

    res = solve(qprob, QuAB3(nreg), dt = h)
    @test isapprox.(s, res, atol = 0.3) |> all

    res = solve(qprob, QuAB4(nreg),dt = h)
    @test isapprox.(s, res, atol = 0.3) |> all
end;
