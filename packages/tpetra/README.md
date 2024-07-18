# Tpetra: Templated Linear Algebra Services Package

Tpetra implements parallel sparse linear algebra and data
redistribution.

"Parallel" includes both MPI (the Message Passing Interface for
distributed-memory parallel processing) and shared-memory threads
(through the Kokkos shared-memory parallel programming model).

"Sparse linear algebra" includes sparse graphs, sparse matrices, block
sparse matrices (where "blocks" are small and dense), and dense
vectors.

"Data redistribution" includes

  - data distributions over MPI processes,
  - reusable communication patterns between data distributions, and
  - hooks to write your own distributed data structures, without
    needing to write MPI communication yourself.

For more details, please refer to the Frequently Asked Questions
document, ./doc/FAQ.txt, or to Tpetra's Doxygen documentation online.


## Copyright and License
See tpetra/COPYRIGHT, tpetra/LICENSE, https://trilinos.github.io/license.html and individual file headers for additional information.


## Questions? 
Contact lead developers:

* Tpetra team     (GitHub handle: @trilinos/tpetra)
