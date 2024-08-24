using Yao
using QuDiffEq
using LinearAlgebra
using OrdinaryDiffEq
using Yao.BitBasis
using Random
using Test
using Yao.YaoBlocks

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

    sol = solve(prob, Tsit5(), dt = 0.1, adaptive = false)
    s = sol.u[end]

    out = solve(qprob, QuLDE(k))

    @test isapprox.(s, out, atol = 0.01) |> all

    qprob = QuLDEProblem(An, b, x, tspan)
    f(u,p,t) = An*u + b;
    prob = ODEProblem(f, x, tspan)

    sol = solve(prob, Tsit5(), dt = 0.1, adaptive = false)
    s = sol.u[end]

    out = solve(qprob, QuLDE(k))
    @test isapprox.(s, out, atol = 0.04) |> all

    # u0 equal to zero
    tspan = (0.0,0.1)

    qprob = QuLDEProblem(Au,b,tspan)
    out = solve(qprob,QuLDE(k))
    t = tspan[2] - tspan[1]
    r_out = (exp(Au*t) - Diagonal(ones(length(b))))*Au^(-1)*b
    @test isapprox.(r_out, out, atol = 1e-3) |> all

    qprob = QuLDEProblem(An,b,tspan)
    out = solve(qprob,QuLDE(k))
    t = tspan[2] - tspan[1]
    r_out = (exp(An*t) - Diagonal(ones(length(b))))*An^(-1)*b
    @test isapprox.(r_out, out, atol = 1e-3) |> all
end
