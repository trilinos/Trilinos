Kokkos implementation of miniFE application version 1.0

In the original miniFE, the only parallel computing that was done was the CG_Solve method. 
Using Kokkos, we have changed this so now the only serial computing is initializing the Mesh structure.

The Kokkos implementation currently only creates a cube of elements.

We have implemented a parallel gather operation to get rid of race conditions.

We gather into a true CRS matrix, allowing problem to be completely general.

Our CG_Solve.hpp is similar to cg_solve.hpp in miniFE 1.0, in that we call functors to allow for 
parallel computing. The only difference is we are using the Kokkos library to run our paralel_for
and parallel_reduce.

Changing the build script will allow you to compile for different devices, similar to minFE 1.0


