# QuDiffEq

[![CI](https://github.com/QuantumBFS/QuDiffEq.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/QuantumBFS/QuDiffEq.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/QuantumBFS/QuDiffEq.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/dgan181/QuDiffEq.jl)

Quantum algorithms for solving differential equations.

This project is part of Julia's Season of Contribution 2019.

For an introduction to the algorithms and an overview of the features, you can take a look at the blog posts: [#1](https://nextjournal.com/dgan181/julia-soc-19-quantum-algorithms-for-differential-equations/edit), [#2](https://nextjournal.com/dgan181/jsoc-19-non-linear-differential-equation-solver-and-simulating-of-the-wave-equation/edit). 

## Installation

<p>
QuDiffEq is a &nbsp;
    <a href="https://julialang.org">
        <img src="https://raw.githubusercontent.com/JuliaLang/julia-logo-graphics/master/images/julia.ico" width="16em">
        Julia Language
    </a>
    &nbsp; package. To install QuDiffEq,
    please <a href="https://docs.julialang.org/en/v1/manual/getting-started/">open
    Julia's interactive session (known as REPL)</a> and press <kbd>]</kbd> key in the REPL to use the package mode, then type the following command
</p>

```julia
pkg> add QuDiffEq
```

## Algorithms

- Quantum Algorithms for Linear Differential Equations,
  - Based on truncated Taylor series
  - Based on HHL.
- Quantum Algorithms for Non Linear Differential Equations.

## Built With

* [Yao](https://github.com/QuantumBFS/Yao.jl) - A framework for Quantum Algorithm Design
* [QuAlgorithmZoo](https://github.com/QuantumBFS/QuAlgorithmZoo.jl) - A repository for Quantum Algorithms

## Authors

See the list of [contributors](https://github.com/QuantumBFS/QuDiffEq.jl/graphs/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/QuantumBFS/QuDiffEq.jl/blob/master/LICENSE) file for details

## References
- D. W. Berry. *High-order quantum algorithm for solving linear differential equations* (https://arxiv.org/abs/1010.2745)
- Tao Xin et al. *A Quantum Algorithm for Solving Linear Differential Equations: Theory and Experiment* (https://arxiv.org/abs/1807.04553)
- Sarah K. Leyton, Tobias J. Osborne. *A quantum algorithm to solve nonlinear differential equations*(https://arxiv.org/abs/0812.4423)
- P. C.S. Costa et al. *Quantum Algorithm for Simulating the Wave Equation* (https://arxiv.org/abs/1711.05394)
