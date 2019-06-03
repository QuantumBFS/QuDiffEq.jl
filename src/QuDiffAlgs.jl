export QuLDE, LDEMSAlgHHL
export QuEuler, QuLeapfrog, QuAB2, QuAB3, QuAB4

abstract type QuODEAlgorithm <: DiffEqBase.AbstractODEAlgorithm end

abstract type LDEMSAlgHHL <: QuODEAlgorithm end

"""
Linear differential equation solvers (non-HHL)

ref : arxiv.org/abs/1807.04553
"""

struct QuLDE <: QuODEAlgorithm end

"""
Linear differential equation solvers using HHL

ref : arxiv.org/abs/1010.2745v2
"""
struct QuEuler{T} <: LDEMSAlgHHL
    step::Int
    α::Vector{T}
    β::Vector{T}

    QuEuler(::Type{T} = Float64) where {T} = new{T}(1,[1.0,],[1.0,])
end

struct QuLeapfrog{T} <: LDEMSAlgHHL
    step::Int
    α::Vector{T}
    β::Vector{T}

    QuLeapfrog(::Type{T} = Float64) where {T} = new{T}(2,[0, 1.0],[2.0, 0])
end
struct QuAB2{T} <: LDEMSAlgHHL
    step::Int
    α::Vector{T}
    β::Vector{T}

    QuAB2(::Type{T} = Float64) where {T} = new{T}(2,[1.0, 0], [1.5, -0.5])
end
struct QuAB3{T} <: LDEMSAlgHHL
    step::Int
    α::Vector{T}
    β::Vector{T}
    QuAB3(::Type{T} = Float64) where {T} = new{T}(3,[1.0, 0, 0], [23/12, -16/12, 5/12])
end
struct QuAB4{T} <: LDEMSAlgHHL
    step::Int
    α::Vector{T}
    β::Vector{T}

    QuAB4(::Type{T} = Float64) where {T} = new{T}(4,[1.0, 0, 0, 0], [55/24, -59/24, 37/24, -9/24])
end
