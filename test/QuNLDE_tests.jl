using Yao
using QuDiffEq
using LinearAlgebra
using BitBasis
using Random
using Test
using YaoBlocks

#Linear Diff Equation Unitary M
function f(du,u,p,t)
    du[1] = -3*u[1]^2 + u[2]
    du[2] = -u[2]^2 - u[1]*u[2]
end
@testset "QuLDE_Test" begin
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

    qprob = QuLDEProblem(A, x, tspan)
    out = transformfunc(qprob.A, qprob.u0, k )
    r_out = zero(x)
    f(r_out, x,1,1)
    @test isapprox.(r_out, 2*out[2:3], atol = 1e-3) |> all
end
