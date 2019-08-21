# Linear Differential Equations

A linear differential equation is written as
```math
 \frac{dx}{dt} = Mx + b
```
 `M` is an arbitrary N × N matrix, `x` and `b` are N dimensional vectors.

`QuDiffEq` allows for the following methods of solving a linear differential equation:
- `QuLDE`: LDE algorithm based on Taylor Truncation. This method evaluates the vector at the last time step, without going through the intermediate steps, unlike other solvers.
- `<: LDEMSAlgHHL`: LDE algorithm based on HHL
  - `QuEuler`
  - `QuLeapfrog`
  - `QuAB2`
  - `QuAB3`
  - `QuAB4`

## Usage

Firstly, we need to define a `QuLDEProblem` for matrix `M` (may be time dependent), initial vector `x` and vector `b` (may be time dependent). `tspan` is the time interval.

```@example lin
using QuDiffEq
using OrdinaryDiffEq, Test
using Random
using LinearAlgebra

siz = 2
M = rand(ComplexF64,siz,siz)
b = normalize!(rand(ComplexF64, siz))
x = normalize!(rand(ComplexF64, siz))
tspan = (0.0,0.4)

qprob = QuLDEProblem(M,b,x,tspan)
```
To solve the problem we use `solve()` after deciding on an algorithm e.g. `alg = QuAB3()` . Here, is an example for `QuLDE`.
```@example lin
alg = QuLDE()
res = solve(qprob,alg)
```
Let's compare the result with a `Tsit5()` from `OrdinaryDiffEq`
```@example lin
f(u,p,t) = M*u + b;
prob = ODEProblem(f, x, tspan)

sol = solve(prob, Tsit5(), dt = 0.1, adaptive = false)
s = sol.u[end]
@test isapprox.(s, res, atol = 0.02) |> all
```
Note : `QuLDE` works only with constant `M` and `b`. There is no such restriction on the other algorithms.
