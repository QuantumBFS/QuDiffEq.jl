# Taylor Truncation

Taylor truncation based Hamiltonian simulation (https://arxiv.org/abs/1412.4687) has many clear advantages. It has better complexity dependence on the precision and allows a greater range of Hamiltonians to be simulated.

The package provides circuits for five kinds of problems:
- Unitary Taylor simulation
- Non-unitary Taylor simulation
- Unitary QuLDE Problem
- Non-unitary QuLDEProblem
- Solution by linearising a non-linear differential equation

To simulate the following the algorithm, simulates the Taylor expansion of ``e^{iHt}`` upto ``k`` orders.
```
|u(t)⟩ = e^{iHt}|u(0)⟩
```
We have,

```math
|x(t)\rangle \approx \sum^{k}_{m=0}\frac{||x(0)||(||M||t)^{m}}{m!}\mathcal{M}^{m}|x(0)\rangle
```
where ``M = iH`` and ``M = ||M||\mathcal{M}``

There are two cases to consider:

1. ``\mathcal{M}`` is unitary. In addition to the vector state register, we have ``T = log_{2}(k+1)`` ancillary bits. The ancilla is in ``|0⟩`` to begin with. The `VS1` block acts on the ancilla register to generate an appropriate superpostion (``\sum^{k}_{m=0}\frac{||x(0)||(||M||t)^{m}}{m!}``). These bits then control the multiplication of ``\mathcal{M}``. ``\mathcal{M}^{j}`` block is control-multiplied by ``|j⟩`` in the ancilla register. `VS1'` is the adjoint of `VS1` and it un-computes the ancilla. After the un-computation we obtain the desired result in all zero ancilla bits.

2. ``\mathcal{M}`` is non-unitary. ``\mathcal{M}`` is expressed as a linear combination of four (at most) unitary  i.e. ``\mathcal{M} =\sum_{i} 1/2 F_{i}``. We have two registers of sizes ``k`` and ``2k``. These registers participate in control-multiplications, as control bits. `VS1` behaves differently to that in the unitary case. In the first register, states with ``j`` ``1``'s (where ``j \in \{0,1,...,k\}``) are raised to probability amplitudes equal to the term with the ``j^{th}`` power in the summation above, while rest of the states are given zero probability. The mappings used is ``m = 2^{k} - 2^{j}`` , ``m`` corresponds to the basis state in the first register. This register governs the the power ``F_i`` need to be raised. The second register is superposed by `VT`, where each new state corresponds to an ``F-i``. When un-computed and measured in the zero state ancilla state, we obtain the desired result.

```@autodocs
Modules = [QuDiffEq]
Pages  = ["TaylorTrunc.jl",]
```

## Quantum Linear Differential Equation

A linear differential equation is written as
```math
 \frac{dx}{dt} = Mx + b
```
 ``M`` is an arbitrary N × N matrix, ``x`` and ``b`` are N dimensional vectors.

The modified version of the `taylorcircuit` is employed here. The transformation over ``x`` is the same as in the `taylorcircuit`, but there is simultaneous transformation over `b` as well. There is an additional ancillary bit that allows the distinction between ``x`` and ``b``. This superpostion is brought about by `V` matrix. Like `VS1` in `taylorcircuit`, `VS2` facilitates the transformation over ``b`` in the simulation.

```@autodocs
Modules = [QuDiffEq]
Pages  = ["QuLDE.jl"]
```

## Quantum Non-linear Differential Equation

The non-linear solver constitutes two sub-routines.
Firstly, the function transform sub-routine, which employs of the `taylorcircuit`. The function transform lets us map ``z`` to ``P(z)``, where ``P`` is a quadratic polynomial.
Secondly, the forward Euler method. The polynomial here is : ``z + hf(z)``. ``f`` is the derivative.

```@autodocs
Modules = [QuDiffEq]
Pages  = ["QuNLDE.jl"]
```
