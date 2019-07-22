using QuDiffEq, Yao
using Test, Random, LinearAlgebra

@testset "HHL_tests" begin
    include("HHL_tests.jl")
end

@testset "PhaseEstimation_tests" begin
    include("PhaseEstimation_tests.jl")
end

@testset "QFT_tests" begin
    include("QFT_tests.jl")
end

@testset "QuLDE_test" begin
    include("QuLDE_tests.jl")
end

@testset "QuNLDE_test" begin
    include("QuNLDE_tests.jl")
end

@testset "QuDiffHHL_tests" begin
    include("QuDiffHHL_tests.jl")
end

@testset "TaylorTrunc_tests" begin
    include("TaylorTrunc_tests.jl")
end
