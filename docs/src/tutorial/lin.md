# Linear Differential Equations

A linear differential equation is written as
```math
 \frac{dx}{dt} = Mx + b
```
 `M` is an arbitrary N × N matrix, `x` and `b` are N dimensional vectors.

`QuDiffEq` allows for the following methods of solving a linear differential equation.
## QuLDE
- `QuLDE`: LDE algorithm based on Taylor Truncation. This method evaluates the vector at the last time step, without going through the intermediate steps, unlike other solvers.

 The exact solution for ``x(t)`` is give by -

```math
 x(t) = e^{Mt}x(0) + (e^{Mt} - I)M^{-1}b
```
This can be Taylor expanded up to the ``k^{th}``order as -
```math
x(t) \approx \sum^{k}_{m=0}\frac{(Mt)^{m}}{m!}x(0) + \sum^{k-1}_{n=1}\frac{(Mt)^{n-1}t}{n!}b
```
The vectors ``x(0)`` and ``b`` are encoded as state - ``|x(0)\rangle = \sum_{i} \frac{x_{i}}{||x||} |i\rangle`` and ``|b\rangle = \sum_{i} \frac{b_{i}}{||b||} |i\rangle`` with ``\{|i \rangle\}``as the computational basis states. We can also write ``M`` as ``M = ||M||\mathcal{M}``. We then get:

```math
|x(t)\rangle \approx \sum^{k}_{m=0}\frac{||x(0)||(||M||t)^{m}}{m!}\mathcal{M}^{m}|x(0)\rangle + \sum^{k-1}_{n=1}\frac{||b||(||M||t)^{n-1}t}{n!}\mathcal{M}^{n-1}|b\rangle
```
To bring about the above transformation, we use the `quldecircuit`.

Note : `QuLDE` works only with constant `M` and `b`. There is no such restriction on the other algorithms.

## LDEMSAlgHHL

 - `QuEuler`
 - `QuLeapfrog`
 - `QuAB2`
 - `QuAB3`
 - `QuAB4`

The HHL algorithm is used for solving a system of linear equations. One can model multistep methods as linear equations, which then can be simulated through HHL.

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
