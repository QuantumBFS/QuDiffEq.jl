using QuDiffEq
using OrdinaryDiffEq, Test
"""
Linearizing non linear ODE and solving it using QuLDE

du/dt = f
u0 -> inital condition for the ode

Jacobian J and vector b is calculated at every time step h.
Δu -> difference for fixed point

Equation input to QuLDE circuit : d(Δu)/dt = J * Δu + b
k -> order of Taylor Expansion in QuLDE circuit
Δu is added to previous value of u.
"""
function f(du,u,p,t)
    du[1] = -2*(u[2]^2)*u[1]
    du[2] = 3*(u[1]^(1.5)) - 0.1*u[2]
end

u0 = [0.2,0.1]
Δu = [1e-6,1e-6]
h = 0.1
k = 3
tspan = (0.0,0.8)
prob = ODEProblem(f,u0,tspan)

qsol = solve(prob,QuLDE(k,Δu),dt = 0.1)
sol = solve(prob,Tsit5(),dt = 0.1,adaptive = false)
v = transpose(hcat(sol.u...))

@test isapprox.(v,qsol, atol = 1e-3) |> all
