using Yao
using QuDiffEq
using LinearAlgebra
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

@testset "TaylorTrunc_Test" begin
    Random.seed!(2)
    N = 1
    k = 3
    tspan = (0.0,0.1)
    Au,An,b,x = diffeqProblem(N)

    qprob = QuLDEProblem(Au, b, x, tspan)
    r,N = taylorsolve(qprob.A,qprob.u0,k,tspan[2])
    out = N*vec(state(r))
    r_out = exp(qprob.A*tspan[2])*qprob.u0
    @test isapprox.(r_out, out, atol = 1e-3) |> all

    qprob = QuLDEProblem(An, b, x, tspan)
    r,N = taylorsolve(qprob.A,qprob.u0,k,tspan[2])
    out = N*vec(state(r))
    r_out = exp(qprob.A*tspan[2])*qprob.u0
    @test isapprox.(r_out, out, atol = 1e-3) |> all
end
