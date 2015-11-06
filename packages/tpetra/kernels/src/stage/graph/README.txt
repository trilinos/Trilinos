To enable compilation see Trilinos/packages/tpetra/kernels/CMakeLists.txt, 
and remove the commented out portion at the bottom.

The package includes couple graph coloring algorithms, and 
a Symmetric Gauss Seidel implementation. The performance tests can be 
found in Trilinos/packages/tpetra/kernels/perf_test/graph 
where there are test for graph coloring, and symetric gauss seidel 
plugging in Preconditioned Conjugate Gradient Algorithm.

