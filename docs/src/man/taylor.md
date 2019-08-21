#Taylor Truncation

Taylor truncation based Hamiltonian simulation (https://arxiv.org/abs/1412.4687) has many clear advantages. It has better complexity dependence on the precision and allows a greater range of Hamiltonians to be simulated.

The package provides circuits for five kinds of problems:
- Unitary Taylor simulation
- Non-unitary Taylor simulation
- Unitary QuLDE Problem
- Non-unitary QuLDEProblem
- Solution by linearising a non-linear differential equation

```@autodocs
Modules = [QuDiffEq]
Pages  = ["TaylorTrunc.jl",]
```

## Quantum linear differential equation solver

A linear differential equation has the form :
```math
 \frac{dx}{dt} = Mx + b
```
``M``is an arbitrary ``N × N`` matrix,``x``and ``b``are dimensional vectors. The exact solution for ``x(t)`` is give by -

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
The `quldecircuit` effects this transformation on the input state to obtain ``x(t)``.

```@autodocs
Modules = [QuDiffEq]
Pages  = ["QuLDE.jl"]
```

## Quantum non-linear differential equation

The non-linear solver constitutes two sub-routines.
Firstly, the function transform sub-routine, which employs of the `taylorcircuit`. The function transform lets us map ``z`` to ``P(z)``, where ``P`` is a quadratic polynomial.
Secondly, the forward Euler method.

```@autodocs
Modules = [QuDiffEq]
Pages  = ["QuNLDE.jl"]
```
