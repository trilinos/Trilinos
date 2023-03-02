The package includes graph coloring algorithms, and an implementation
of symmetric multithreaded Gauss-Seidel that uses those graph coloring
algorithms. 

Both algorithms have performance tests in
Trilinos/packages/tpetra/kernels/perf_test/graph.  The performance
test for symmetric Gauss-Seidel exercises the algorithm in the context
of a preconditioned iterative linear solve, using the Method of
Conjugate Gradients (CG).

