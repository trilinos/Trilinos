# ShyLU: Scalable Hybrid LU Preconditioner and Solver

ShyLU is a package for solving sparse linear systems using domain decomposition methods. ShyLU has two main focus areas - 
(1) distributed memory domain-decomposition solvers and (2) node-level solvers and kernels that support the distributed memory solvers. The approaches in ShyLU are algebraic and so can be used as a black-box solvers.

The node level solvers include sparse LU factorization (basker), sparse Cholesky factorization (Tacho), multithreaded triangular solver (HTS) and a fast iterative ILU algorithm (FastILU).

## Domain Decomposition Solvers

The Domain Decomposition methods in ShyLU are a one-level hybrid direct/iterative approach based on Schur complements, two-level balancing domain decomposition method (bddc) and Overlapping Schwarz methods (FROSch). 

The one-level hybrid direct/iterative approach based on Schur complements was developed to provide robustness similar to sparse direct solvers with less memory use. This has been particularly effective for circuit simulation problems.

[FROSch](https://shylu-frosch.github.io) can be used for both one-level overlapping Schwarz and two-level GDSW (Generalized Dryja-Smith-Widlund) type preconditioners. FROSch has been shown to be effective in several applications including land-ice simulations.

BDDC is an experimental code for the solver that implements the Balancing Domain Decomposition solver.

## Node Level Direct Solvers

Basker is a shared-memory parallel LU factorization based direct solver that uses the BTF ordering.

Tacho is a shared-memory parallel Cholesky factorization based direct solver that uses Kokkos based tasking.

HTS is a sparse triangular solver that uses OpenMP for triangular solves on the host.

FastILU is an iterative ILU and triangular solve implementation using Kokkos.

## Documentation

ShyLU is part of the [Trilinos Project](https://trilinos.github.io), and additional information (e.g., examples, tutorials, and source code documentation) is available through [ShyLU's Doxygen webpages](https://trilinos.github.io/docs/shylu/index.html).

## Questions? 

Contact lead developers:

* **Shylu team**:          GitHub handle: @trilinos/shylu
* **Ichitaro Yamazaki**:   GitHub handle: [@iyamazaki](https://github.com/iyamazaki) or iyamaza@sandia.gov
* **Nathan Ellingwood**:   GitHub handle: [ndellingwood](https://github.com/ndellingwood)
* **Siva Rajamanickam**:   GitHub handle: [srajama1](https://github.com/srajama1) or srajama@sandia.gov

## Copyright and License

For general copyright and license information, refer to the Trilinos [License and Copyright](https://trilinos.github.io/about.html#license-and-copyright) page.

For ShyLU-specific copyright and license details, refer to the [shylu/COPYRIGHT](COPYRIGHT) and [shylu/LICENSE](LICENSE) files located in the `shylu` directory. Additional copyright information may also be found in the headers of individual source files.

For developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.
