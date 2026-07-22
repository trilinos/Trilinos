# Adelus

## Overall Description

Adelus performs LU factorization with partial pivoting and solves
for a dense (real or complex) linear equation system on a distributed
computing system using MPI for message passing. It can be considered
a replacement for Pliris, which runs only on CPUs. Adelus uses
Kokkos and Kokkos Kernels to achieve performance portability on
heterogeneous architectures equipped with CPUs and accelerated
hardware such as GPUs.

The matrix is blocked mapped onto MPI processes (either on host
memory or device memory). It is factored and solved as if it was
torus-wrapped for an optimal balance of computation and communication.
A permutation operation is performed on the results to restore the
results to their correct position so the torus-wrap distribution
is hidden to the user. Each process contains blocks of the matrix
and the right hand sides. A distribution function is used to optimally
load balance the computation and communication during the LU
factorization of the matrix. Each process handles its own workload
through the Kokkos Kernels BLAS routines optimized for multi-threaded
CPUs and massively parallel GPUs. GPU-aware MPI can be used on GPU
architectures to allows direct communication among GPU buffers.

Adelus provides three interfaces, FactorSolve (factors and solves
a dense system in which matrix and RHS are packed in Kokkos View),
Factor (factors a dense matrix for later solve), and Solve (solves
a previously factored dense matrix for provided RHS). Adelus supports
solving for multiple RHS vectors. Adelus allows applications to
concurrently launch multiple solvers each of which is associated
with a MPI sub-communicator of the MPI_COMM_WORLD.

### Relation of the Matrix Distribution to MPI Process Topology

The matrix is distributed to the different MPI processes using the following criterion. 
No MPI process can have no more than 1 row or column than any other process. This in effect 
balances the the number of matrix entries that each process needs to fill. An example 
illustrating this algorithm follows.

```
 Example: 6 MPI processes
        nprocsr =3  (Number of processes for a row)
        Process id in the box below.


          0      1      2      < --- my_col
       ------ ------ ------
       |    | |    | |    |
       |  0 | |  1 | |  2 |    < --- my_row = 0 for these processes
       |    | |    | |    |
       ------ ------ ------
       |    | |    | |    |
       |  3 | |  4 | |  5 |    < --- my_row = 1 for these processes
       |    | |    | |    |
       ------ ------ ------

    For an 1000 x 1000 matrix

      Process id           0        1       2      3     4     5
      my_rows             334      333     333    334   333   333
      my_cols             500      500     500    500   500   500
      my_first_row         1       335     668     1    335   668
      my_first_col         1        1       1     501   501   501

      Note: The fortran convention is assumed the matrix begins with index 1.
```
In addition, the Right Hand Sides (RHS) are distributed in a similar manner, i.e.
no MPI process can have no more than 1 RHS than any other process. Using the example 
shown above the RHS are distributed:

   For 1 RHS (x)

         0         1      2        < --- my_col
       ------    ------ ------
       |    |x   |    | |    |
       |  0 |x   |  1 | |  2 |     < --- my_row = 0 for these processes
       |    |x   |    | |    |
       ------    ------ ------
       |    |x   |    | |    |
       |  3 |x   |  4 | |  5 |     < --- my_row = 1 for these processes
       |    |x   |    | |    |
       ------    ------ ------

   For 4 RHS (x)

         0            1        2        < --- my_col
       ------     ------     ------
       |    |xx    |    |x   |    |x
       |  0 |xx    |  1 |x   |  2 |x     < --- my_row = 0 for these processes
       |    |xx    |    |x   |    |x
       ------     ------     ------
       |    |xx    |    |x   |    |x
       |  3 |xx    |  4 |x   |  5 |x     < --- my_row = 1 for these processes
       |    |xx    |    |x   |    |x
       ------     ------     ------

This first version of the solver requires the RHS before the solve occurs since
the forward solve is realized by factoring the matrix with the RHS appended to 
the matrix.


### Directory Organization

We organize the directories as follows:

1. Public interfaces to the solver and the distribution function live in the 
```src/``` subdirectory (```adelus/src```):

* ```Adelus::GetDistribution()```: gives the distribution information that is required
by the dense solver to the user that defines the matrix block and right hand side information.

* ```Adelus::AdelusHandle<...>```: an application must create a handle to the Adelus communicator and necessary metadata (the handle is passed to every subsequent Adelus function call)

* ```Adelus::FactorSolve()```: factors and solves the dense matrix in which the matrix
and rhs are packed in Kokkos View

* ```Adelus::FactorSolve_devPtr()```: matrix and rhs are packed and passed as device pointer

* ```Adelus::FactorSolve_hostPtr()```: matrix and rhs are packed and passed as host pointer

* ```Adelus::Factor()```: factors the dense matrix for later solve

* ```Adelus::Solve()```: solves the previously factored dense matrix for provided RHS

2. Implementations of the phases of the solver (i.e. factor, solve, permutation)  
and other utility functions also locate in the ```src/``` subdirectory.

3. Correctness tests is in the ```test/``` subdirectory.

4. A simple example that generates a random matrix and a right-hand-side to
    exercise the solver is in the ```example/``` subdirectory.


## Configuring, Building, and Installing Adelus

 Adelus mostly consists of header files. Only a few functions contained in .cpp files
have to be compiled into object files outside of the application's source code. It
should be noted that a C++11 compliant compiler is needed to build Adelus. Since
Adelus is distributed within Trilinos, and uses Kokkos and Kokkos Kernels extensively,
Trilinos' CMake build system is preferred. To enable Adelus when building Trilinos,
set the CMake option ```Trilinos_ENABLE_Adelus```. Trilinos' build system lets packages
express dependencies on other packages or external libraries. Adelus has a required
dependency on Kokkos and Kokkos Kernels, Trilinos will enable Kokkos and Kokkos Kernels
automatically. Following the Kokkos and Kokkos Kernels style, Adelus's features are
controlled via CMake options in the form of ```Adelus_ENABLE_OPTION```. A list of Adelus
options can be found below.

* ```Adelus_ENABLE_ZCPLX```
  * Whether to enable double precision complex functionality
  * ```BOOL``` Default: ON
* ```Adelus_ENABLE_SCPLX```
  * Whether to enable single precision complex functionality
  * ```BOOL``` Default: OFF
* ```Adelus_ENABLE_DREAL```
  * Whether to enable double precision functionality
  * ```BOOL``` Default: OFF
* ```Adelus_ENABLE_SREAL```
  * Whether to enable single precision functionality
  * ```BOOL``` Default: OFF
* ```Adelus_ENABLE_TIMING```
  * Whether to enable internal solver timing
  * ```BOOL``` Default: OFF
* ```Adelus_ENABLE_HOSTPINNEDMEM```
  * Whether to use Cuda/HIP Host Pinned memory for MPI
  * ```BOOL``` Default: OFF
* ```Adelus_ENABLE_USEDEEPCOPY```
  * Whether to Use Kokkos::deep_copy for BLAS copy
  * ```BOOL``` Default: OFF
* ```Adelus_ENABLE_PRINTSTATUS```
  * Whether to enable status prints
  * ```BOOL``` Default: OFF

 We refer readers to Trilinos', Kokkos', and Kokkos Kernels' documentations for
further details of building Trilinos, Kokkos, and Kokkos Kernels.

 Below is an example of a Trilinos build script enabling Adelus (and Kokkos and
Kokkos Kernels) with CUDA back-end (NVIDIA Volta generation CC 7.0) and OpenMP host
back-end to perform comparison between the CUDA code and the more classical MPI
version of the code.  The example below includes the environment variables and Cmake
options to build Adelus.

### Configure

```
export TRILINOS_HOME=/your/Trilinos/home/directory
export CXX=${TRILINOS_HOME}/Trilinos/packages/kokkos/bin/nvcc_wrapper
export CUDA_LAUNCH_BLOCKING=1
export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
export OMPI_CXX=${TRILINOS_HOME}/Trilinos/packages/kokkos/bin/nvcc_wrapper
export LLNL_USE_OMPI_VARS="Y"
export NVCC_WRAPPER_DEFAULT_COMPILER=`which g++`

INSTALL_LOCATION=${TRILINOS_HOME}/install-trilinos

cmake \
\
-D CMAKE_INSTALL_PREFIX:PATH="${INSTALL_LOCATION}" \
-D CMAKE_CXX_COMPILER:FILEPATH="`which ${CXX}`" \
-D CMAKE_C_COMPILER:FILEPATH="`which mpicc`" \
-D CMAKE_Fortran_COMPILER:FILEPATH="`which mpif77`" \
-D CMAKE_CXX_FLAGS="-O3 -g -Wall -Wno-unknown-pragmas -Wno-unused-but-set-variable -Wno-inline -Wshadow -I${MPI_ROOT}/include" \
-D CMAKE_C_FLAGS="-O3 -g" \
-D CMAKE_Fortran_FLAGS="-g" \
\
-D TPL_ENABLE_MPI:BOOL=ON \
-D TPL_ENABLE_Matio:BOOL=OFF \
-D TPL_ENABLE_X11:BOOL=OFF \
-D TPL_ENABLE_BLAS:BOOL=ON \
-D TPL_ENABLE_CUDA:BOOL=ON \
-D BLAS_INCLUDE_DIRS="${BLAS_HEADER_DIR}" \
-D BLAS_LIBRARY_DIRS="${BLAS_LIB_DIR}" \
-D BLAS_LIBRARY_NAMES="blas_lib_name" \
\
-D MPI_BASE_DIR="${MPI_ROOT}" \
\
-D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
-D Trilinos_ENABLE_CXX11:BOOL=ON \
-D Trilinos_ENABLE_Kokkos:BOOL=ON \
-D Trilinos_ENABLE_Adelus:BOOL=ON \
-D Trilinos_ENABLE_TESTS:BOOL=ON \
-D Trilinos_ENABLE_EXAMPLES:BOOL=OFF \
-D Trilinos_ENABLE_COMPLEX_DOUBLE:BOOL=ON \
-D Trilinos_EXTRA_LINK_FLAGS:STRING="-lmpi_ibm -fopenmp" \
\
-D Kokkos_ENABLE_SERIAL:BOOL=OFF \
-D Kokkos_ENABLE_OPENMP:BOOL=ON \
-D Kokkos_ENABLE_THREADS:BOOL=OFF \
-D Kokkos_ENABLE_CUDA:BOOL=ON \
-D Kokkos_ENABLE_CUDA_UVM:BOOL=OFF \
-D Kokkos_ENABLE_CUDA_LAMBDA:BOOL=ON \
-D Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE:BOOL=ON \
-D Kokkos_ENABLE_TESTS:BOOL=OFF \
-D Kokkos_ARCH_POWER9=ON \
-D Kokkos_ARCH_VOLTA70=ON \
\
-D Adelus_ENABLE_ZCPLX:BOOL=ON \
-D Adelus_ENABLE_TIMING:BOOL=ON \
-D Adelus_ENABLE_HOSTPINNEDMEM:BOOL=OFF \
-D Adelus_ENABLE_PRINTSTATUS:BOOL=ON \
\
${TRILINOS_HOME}/Trilinos
```

### Build

```
make -j N 
```

### Install

```
make -j N install
```

## Using Adelus Test Driver

 This is an example that drives Adelus to solve for a dense linear system with
one RHS vector:

1. Determine the matrix and RHS vector distribution on each MPI process using
the ```Adelus::GetDistribution``` utility function.

2. Once the portions of matrix and RHS vectors are computed on each process,
the solver can be called. In this example, the portion of matrix on each MPI
process and the reference solution vector are randomly generated. Then, the
assigned RHS vectors on MPI processes can be computed.

3. Create a handle to the Adelus communicator and necessary metadata

4. Launch Adelus using ```Adelus::FactorSolve```, or ```Adelus::FactorSolve_devPtr```,
or ```Adelus::FactorSolve_hostPtr```.

5. Gather results.

6. Compare the returned solution vector with the reference vector.

### Compile with Makefile

```
make -j
```

### Compile with CMake

Make a ```cmakebuild``` directory, and then run the command:

```
cmake -DTrilinos_DIR={trilinos_install_path}/include ..
```

from the ```cmakebuild``` directory.

## Documentation

Additional information on Adelus can be found at the [Trilinos Project](https://trilinos.github.io), and through the [Adelus's Doxygen webpages](https://trilinos.github.io/docs/adelus/index.html).

## Questions?

Any questions can be directed to: Vinh Dang (vqdang@sandia.gov), Joseph Kotulski (jdkotul@sandia.gov), Siva Rajamanickam (srajama@sandia.gov)

## Copyright and License

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Adelus is licensed under standard 3-clause BSD terms of use.
For Adelus-specific copyright and license details, refer to the [adelus/COPYRIGHT](COPYRIGHT) and [adelus/LICENSE](LICENSE) files located in the `adelus` directory. Additional copyright information may also be found in the headers of individual source files.

For general copyright and license information, refer to the Trilinos [License](../../LICENSE) and [Copyright](../../COPYRIGHT).

For developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.
