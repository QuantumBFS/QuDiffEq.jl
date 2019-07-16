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
h = 0.1
k = 2
tspan = (0.0,0.8)
prob = ODEProblem(f,u0,tspan)

qsol = solve(prob,QuLDE(k),dt = h)
sol = solve(prob,Tsit5(),dt = h,adaptive = false)
r_out = transpose(hcat(sol.u...))

@test isapprox.(r_out,qsol, atol = 1e-3) |> all
