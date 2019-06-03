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
    A = rand_unitary(siz)
    b = normalize!(rand(ComplexF64, siz))
    x = normalize!(rand(ComplexF64, siz))
    A, b, x
end

@testset "QuLDE_Test" begin
    Random.seed!(2)
    N = 1
    k = 3
    tspan = (0.0,0.4)
    A,b,x = diffeqProblem(N)

    qprob = QuLDEProblem(A, b, x, tspan)
    f(u,p,t) = A*u + b;
    prob = ODEProblem(f, x, tspan)

    sol = solve(prob, Tsit5(), dt = 0.1, adaptive = :false)
    s = vcat(sol.u...)

    out = solve(qprob, QuLDE(), k)

    @test isapprox.(s[end-1:end], out, atol = 0.01) |> all
end
