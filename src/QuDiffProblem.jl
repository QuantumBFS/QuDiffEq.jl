export QuLDEProblem
#abstract type QuODEProblem{uType,tType,isinplace} <: DiffEqBase.AbstractODEProblem{uType,tType,isinplace} end
struct QuLDEProblem{F,C,U,T} #<: QuODEProblem{uType,tType,isinplace}
    A::F
    b::C
    u0::U
    tspan::NTuple{2,T}

    #function QuLDEMSProblem(A,b,u0,tspan)
    #  new{typeof(u0),typeof(tspan),false,typeof(A),typeof(b)}(A,b,u0,tspan)
    #end
end
