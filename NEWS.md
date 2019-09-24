[QuDiffEq.jl](https://github.com/QuantumBFS/QuDiffEq.jl) was developed as part of JSoC 2019. It allows one to solve differential equations using quantum algorithms.
We have the following algorithms and example in place:
* `QuLDE` for linear differential equations. There is also a rountine for using `QuLDE` to solve non-linear differential equations by making linear approximations of the functions at hand. [Read more.](https://nextjournal.com/dgan181/julia-soc-19-quantum-algorithms-for-differential-equations/edit)
* Several HHL-based multistep methods for linear differential equations.
* `QuNLDE` for solving non-linear quadratic differential equations. [Read more.](https://nextjournal.com/dgan181/jsoc-19-non-linear-differential-equation-solver-and-simulating-of-the-wave-equation/edit)
* Simulation of the wave equation using quantum algorithms. [Read more.](https://nextjournal.com/dgan181/jsoc-19-non-linear-differential-equation-solver-and-simulating-of-the-wave-equation/edit)
