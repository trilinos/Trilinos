![](doc/images/tpetra.jpg)

# Tpetra: Parallel sparse linear algebra

Trilinos’ Tpetra package implements _linear algebra objects_. These include sparse graphs, sparse matrices, and dense vectors. Many Trilinos packages and applications produce, modify, or consume Tpetra’s linear algebra objects, or depend on Tpetra’s parallel data redistribution facilities.

## Parallel

Tpetra is “hybrid parallel,” meaning that it uses at least two levels of parallelism:

*   MPI (the Message Passing Interface) for distributed-memory parallelism
*   Any of various shared-memory parallel programming models within an MPI process

Tpetra uses the Kokkos package’s shared-memory parallel programming model for data structures and computational kernels. Kokkos makes it easy to port Tpetra to new computer architectures, and to extend its use of parallel computational kernels and thread-scalable data structures. Kokkos supports several different parallel programming models, including

*   OpenMP
*   POSIX Threads (Pthreads)
*   Nvidia’s CUDA programming model for graphics processing units (GPUs)

## Distributed

Tpetra also understands parallel distributions of data. It has objects called Maps that let users represent a particular distribution of data over one or more parallel processes. Users may also construct data migrations (Import or Export) between two distributions, and use them either on Tpetra native objects, or on their own objects that implement a simple data packing and unpacking interface (DistObject).

## Features

Tpetra has the following unique features:

*   Native support for representing and solving very large graphs, matrices, and vectors. “Very large” means over two billion unknowns or other entities.
*   Matrices and vectors may contain many different kinds of data, such as floating-point types of different precision, complex-valued types, automatic differentiation objects from the Sacado package, or stochastic PDE discretization types from the Stokhos package.
*   Support for many different shared-memory parallel programming models

If you are familiar with the Epetra package, think of Tpetra as the next version of Epetra. Epetra’s interface is stable, but Tpetra has better support for shared-memory parallelism. Current and future software development effort will focus on Tpetra and the solvers and preconditioners that depend on it.
For more details, please refer to the Frequently Asked Questions
document, [FAQ](doc/FAQ.txt), or to Tpetra's Doxygen documentation online.

## Questions?

Contact the lead developers:

- **Tpetra team**:   GitHub handle: @trilinos/tpetra
- **Chris Siefert**: GitHub handle [csiefer2](https://github.com/csiefer2), email: csiefer@sandia.gov

## Copyright and License

For general copyright and license information, refer to the Trilinos [License and Copyright](https://trilinos.github.io/about.html#license-and-copyright) page.

For Tpetra-specific copyright and license details, refer to the [tpetra/COPYRIGHT](COPYRIGHT) and [tpetra/LICENSE](LICENSE) files located in the `tpetra` directory. Additional copyright information may also be found in the headers of individual source files.

For developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.
