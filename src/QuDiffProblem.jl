export QuLDEProblem
"""
Linear ODE Problem definition

"""
struct QuODEFunction{iip,MType}<: DiffEqBase.AbstractODEFunction{iip}
    linmatrix::Array{MType,2}
    QuODEFunction(linmatrix::Array{MType,2},::Type{PType} = Float32) where {MType,PType} = new{true,Complex{PType}}(linmatrix)
end

abstract type QuODEProblem{uType,tType,isinplace} <: DiffEqBase.AbstractODEProblem{uType,tType,isinplace} end

struct QuLDEProblem{uType,tType,isinplace, F, P} <: QuODEProblem{uType,tType,isinplace}
    A::F
    b::P
    u0::uType
    tspan::Tuple{tType,tType}
    function QuLDEProblem(A::QuODEFunction{iip,MType},b::Array{T,1},u0::Array{G,1},tspan;kwargs...) where {iip,MType,T,G}
        PType = eltype(A.linmatrix)
     new{Array{Complex{PType},1},typeof(tspan),iip,typeof(A.linmatrix),Array{Complex{PType},1}}(A.linmatrix,b,u0,tspan)
    end
    function QuLDEProblem(A,b::Array{T,1},u0::Array{G,1},tspan,::Type{PType} = Float32;kwargs...,) where {T,G,PType}
        f = QuODEFunction(A,PType)
     new{Array{Complex{PType},1},typeof(tspan[1]),isinplace(f),typeof(f.linmatrix),Array{Complex{PType},1}}(f.linmatrix,b,u0,tspan)
    end
end
