using Yao
using QuDiffEq
using LinearAlgebra
using OrdinaryDiffEq
using BitBasis
using Random
using Test
using YaoBlocks

#Linear Diff Equation Unitary M
function diffeqProblem(nbit::Int)
    siz = 1<<nbit
    Au = rand_unitary(siz)
    An = rand(ComplexF64,siz,siz)
    b = normalize!(rand(ComplexF64, siz))
    x = normalize!(rand(ComplexF64, siz))
    Au,An,b, x
end

@testset "QuLDE_Test" begin
    Random.seed!(2)
    N = 1
    k = 3
    tspan = (0.0,0.4)
    Au,An,b,x = diffeqProblem(N)
    qprob = QuLDEProblem(Au, b, x, tspan)
    f(u,p,t) = Au*u + b;
    prob = ODEProblem(f, x, tspan)

    sol = solve(prob, Tsit5(), dt = 0.1, adaptive = :false)
    s = vcat(sol.u[end])

    out = solve(qprob, QuLDE(), k)

    @test isapprox.(s, out, atol = 0.01) |> all

    qprob = QuLDEProblem(An, b, x, tspan)
    f(u,p,t) = An*u + b;
    prob = ODEProblem(f, x, tspan)

    sol = solve(prob, Tsit5(), dt = 0.1, adaptive = :false)
    s = vcat(sol.u[end])

    out = solve(qprob, QuLDE(), k)
    @test isapprox.(s, out, atol = 0.02) |> all
end
