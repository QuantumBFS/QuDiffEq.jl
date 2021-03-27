module QuDiffEq

using Yao
using YaoBlocks
using DiffEqBase
using BitBasis
using LinearAlgebra
using ForwardDiff
using YaoExtensions

include("QuDiffProblem.jl")
include("QuDiffAlgs.jl")
include("TaylorTrunc.jl")
include("QuDiffHHL.jl")
include("HHL.jl")
include("PhaseEstimation.jl")
include("QuLDE.jl")
include("QuNLDE.jl")
end # module
