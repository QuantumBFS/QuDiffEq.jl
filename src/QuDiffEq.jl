module QuDiffEq

using Yao
using Yao.YaoBlocks
using Yao.YaoArrayRegister
using Yao.BitBasis
using DiffEqBase
using LinearAlgebra
using ForwardDiff
import Yao.YaoArrayRegister: u1rows!

include("QuDiffProblem.jl")
include("QuDiffAlgs.jl")
include("TaylorTrunc.jl")
include("QuDiffHHL.jl")
include("HHL.jl")
include("PhaseEstimation.jl")
include("QuLDE.jl")
include("QuNLDE.jl")

end # module
