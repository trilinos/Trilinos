# Belos: Block Linear Solvers Package

## Overview

Belos provides next-generation iterative linear solvers and a powerful linear solver developer framework. This framework includes the following abstract interfaces and implementations:

*   Abstract interfaces to linear algebra using traits mechanisms. This allows the user to leverage any existing investment in their description of matrices and vectors. The provided concrete linear algebra adapters enable Belos to be used anywhere Epetra and Thyra are employed for linear algebra services.
*   Abstract interfaces to orthogonalization; implementations of iterated classical Gram-Schmidt (ICGS), classical Gram-Schmidt with a DGKS correction step, and iterated modified Gram-Schmidt (IMGS) are included.
*   Abstract interfaces to iteration kernels; implementations of conjugate gradient (CG), block CG, block GMRES, pseudo-block GMRES, block flexible GMRES, and GCRO-DR iterations are included.
*   Powerful solver managers are provided for solving a linear system using CG or block CG, GMRES or block GMRES with restarting, pseudo-block GMRES for performing single-vector GMRES simultaneously on multiple right-hand sides, and a single-vector recycled Krylov method (GCRO-DR).
*   Basic linear problem class is provided for the user to define a unpreconditioned or preconditioned (left, right, two-sided) linear system for Belos to solve.

## Publications

If you use Belos in your applications, please cite Belos using the following publication:

*   Amesos2 and Belos: Direct and iterative solvers for large sparse linear systems, Eric Bavier, Mark Hoemmen, Sivasankaran Rajamanickam, Heidi Thornquist, Scientific Programming, 2012, Volume 20, Issue 3.

### Other Publications

*   “[Amesos2 and Belos: Direct and iterative solvers for large sparse linear systems.](https://onlinelibrary.wiley.com/doi/10.3233/SPR-2012-0352)” Eric Bavier, Mark Hoemmen, Sivasankaran Rajamanickam, and Heidi Thornquist. Scientific Programming, 2012.
*   [A Communication-Avoiding, Hybrid-Parallel, Rank-Revealing Orthogonalization Method](https://ieeexplore.ieee.org/document/6012905) Mark Hoemmen. IEEE International Parallel and Distributed Processing Symposium, May 2011.
*   [Cooperative Application/OS DRAM fault recovery](http://dx.doi.org/10.1007/978-3-642-29740-3_28) Patrick G. Bridges, Michael A. Heroux, Mark Hoemmen, Kurt Ferreira, Philip Soltero, and Ronald B. Brightwell. Workshop on Resiliency in High-Performance Computing (Resilience 2011) in conjunction with the 17th International European Conference on Parallel and Distributed Computing (Euro-Par 2011), Bordeaux, France, 29 August -- 02 September 2011.

## Documentation

Belos is part of the [Trilinos Project](https://trilinos.github.io), and additional information (e.g., examples, tutorials, and source code documentation) is available through [Belos's Doxygen webpages](https://trilinos.github.io/docs/belos/index.html).

## Questions?

Contact the lead developers:

- **Belos team**     (GitHub handle: @trilinos/belos)
- **Heidi Thornquist** (GitHub handle: [hkthorn](https://github.com/hkthorn) or hkthorn@sandia.gov)

## Copyright and License

For general copyright and license information, refer to the Trilinos [License and Copyright](https://trilinos.github.io/about.html#license-and-copyright) page.

For Belos-specific copyright and license details, refer to the [belos/COPYRIGHT](COPYRIGHT) and [belos/LICENSE](LICENSE) files located in the `belos` directory. Additional copyright information may also be found in the headers of individual source files.

For developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.

## Copyright and License
See belos/COPYRIGHT, belos/LICENSE, https://trilinos.github.io/license.html and individual file headers for additional information.


## Questions? 
Contact lead developers:

