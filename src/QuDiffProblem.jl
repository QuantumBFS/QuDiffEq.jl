export QuLDEProblem
"""
Linear ODE Problem definition

"""
struct QuODEFunction{iip,MType}<: DiffEqBase.AbstractODEFunction{iip}
    linmatrix::Array{MType,2}
    QuODEFunction(linmatrix::Array{MType,2}) where MType = new{true,ComplexF64}(linmatrix)
end

abstract type QuODEProblem{uType,tType,isinplace} <: DiffEqBase.AbstractODEProblem{uType,tType,isinplace} end

struct QuLDEProblem{uType,tType,isinplace, F, P} <: QuODEProblem{uType,tType,isinplace}
    A::F
    b::P
    u0::uType
    tspan::Tuple{tType,tType}

    function QuLDEProblem(A::QuODEFunction{iip,MType},b::Array{ComplexF64,1},u0::Array{ComplexF64,1},tspan;kwargs...) where {iip,MType}
     new{typeof(u0),typeof(tspan),iip,typeof(A),typeof(b)}(A,b,u0,tspan)
    end
    function QuLDEProblem(A,b::Array{T,1},u0::Array{G,1},tspan;kwargs...) where {T,G}
        f = QuODEFunction(A)
     new{Array{ComplexF64,1},typeof(tspan[1]),isinplace(f),typeof(A),Array{ComplexF64,1}}(f.linmatrix,b,u0,tspan)
    end
end
