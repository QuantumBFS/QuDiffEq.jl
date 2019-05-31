using QuDiffEq, Yao
using Test, Random, LinearAlgebra

#@testset "HHL_tests" begin
#    include("HHL_tests.jl")
#end

@testset "QuLDE_test" begin
    include("QuLDE_tests.jl")
end

@testset "QuDiffHHL_tests" begin
    include("QuLDE_tests.jl")
end
