#  Non-linear Differential Equations

`QuDiffEq` has two algorithms for solving non-linear differential equations:
- `QuNLDE`: Uses function transformation and the forward Euler method.
- `QuLDE`: Linearises the differential equation at every iteration.

## Usage
Let's say we want to solve the following set of differential equations.
```math
\begin{array}{rcl} \frac{dz_1}{dt} & = & z_2 - 3 z_{1}^{2} \\ \frac{dz_2}{dt}  &=& -z_{2}^{2} - z_1 z_{2} \end{array}
```
Let's take the time interval to be from 0.0 to 0.4. We define the in initial vector randomly.
```julia
using QuDiffEq
using OrdinaryDiffEq
using Random
using LinearAlgebra

tspan = (0.0,0.4)
x = [0.6, 0.8]
```

- For `QuNLDE`, we need to define a `<: QuODEProblem`. At present, we use only `QuLDEProblem` as a Qu problem wrapper.
  `QuNLDE` can solve only quadratic differential equations. `A` is the coefficient matrix for the quadratic differential equation.

```julia
N = 2 # size of the input vector
siz = nextpow(2, N + 1)

A = zeros(ComplexF32,2^(siz),2^(siz));
A[1,1] = ComplexF32(1);
A[5,3] = ComplexF32(1);
A[5,6] = ComplexF32(-3);
A[9,11] = ComplexF32(-1);
A[9,7] = ComplexF32(-1);
```

```julia
qprob = QuLDEProblem(A,x,tspan);

```
To solve the problem we use `solve()`
```julia
res = solve(qprob,QuNLDE(), dt = 0.1);
```
Comparing the result with `Euler()`

```julia eval = false

function f(du,u,p,t)
    du[1] = -3*u[1]^2 + u[2]
    du[2] = -u[2]^2 - u[1]*u[2]
end

prob = ODEProblem(f, x, tspan)
sol = solve(prob, Euler(), dt = 0.1, adaptive = false)

using Plots;

plot(sol.t,real.(res),lw = 1,label="QuNLDE()")
plot!(sol,lw = 3, ls=:dash,label="Euler()")
```

![](docs/assets/figures/QuNLDE-plot.svg)

- For `QuLDE`, the problem is defined as a `ODEProblem`, similar to that in OrdinaryDiffEq.jl . `f` is the differential equation written symbolically. We can use prob from the previous case itself.

```julia
res = solve(prob,QuLDE(),dt = 0.1)
```
```
sol = solve(prob, Tsit5(), dt = 0.1, adaptive = false)

using Plots
plot(sol.t,real.(res),lw = 1,label="QuNLDE()")
plot!(sol,lw = 3, ls=:dash,label="Tsit5()")
```
![](docs/assets/figures/QuLDE-plot.svg)
