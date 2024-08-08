[![OpenSSF Scorecard](https://api.securityscorecards.dev/projects/github.com/kokkos/kokkos-kernels/badge)](https://securityscorecards.dev/viewer/?uri=github.com/kokkos/kokkos-kernels)
[![Generic badge](https://readthedocs.org/projects/kokkos-kernels/badge/?version=latest)](https://kokkos-kernels.readthedocs.io/en/latest/)

![KokkosKernels](https://avatars2.githubusercontent.com/u/10199860?s=200&v=4)

# Kokkos Kernels

Kokkos C++ Performance Portability Programming EcoSystem: Math Kernels -
Provides BLAS, Sparse BLAS and Graph Kernels 

KokkosKernels implements local computational kernels for linear
algebra and graph operations, using the Kokkos shared-memory parallel
programming model.  "Local" means not using MPI, or running within a
single MPI process without knowing about MPI.  "Computational kernels"
are coarse-grained operations; they take a lot of work and make sense
to parallelize inside using Kokkos.  KokkosKernels can be the building
block of a parallel linear algebra library like Tpetra that uses MPI
and threads for parallelism, or it can be used stand-alone in your
application.

Computational kernels in this subpackage include the following:

* (Multi)vector dot products, norms, and updates (AXPY-like
    operations that add vectors together entry-wise)
* Sparse matrix-vector multiply and other sparse matrix / dense
    vector kernels
* Sparse matrix-matrix multiply
* Graph coloring
* Gauss-Seidel with coloring (generalization of red-black)
* Other linear algebra and graph operations

We organize this directory as follows:

1. Public interfaces to computational kernels live in the src/
     subdirectory (kokkos-kernels/src):

*    Kokkos_Blas1_MV.hpp: (Multi)vector operations that
       Tpetra::MultiVector uses
*    Kokkos_Sparse_CrsMatrix.hpp: Declaration and definition of
       KokkosSparse::CrsMatrix, the sparse matrix data structure used
       for the computational kernels below
*    KokkosSparse_spmv.hpp: Sparse matrix-vector multiply with a
       single vector, stored in a 1-D View + Sparse matrix-vector multiply with
       multiple vectors at a time (multivectors), stored in a 2-D View

2. Implementations of computational kernels live in the src/impl/
     subdirectory (kokkos-kernels/src/impl)

3. Correctness tests live in the unit_test/ subdirectory, and
     performance tests live in the perf_test/ subdirectory

4. Simple example scripts to build Kokkoskernels are in
     example/buildlib/


Do NOT use or rely on anything in the `KokkosBlas::Impl` namespace, or
on anything in the impl/ subdirectory.

This separation of interface and implementation lets the interface
assign the users' Views to View types with the desired attributes
(e.g., read-only, RandomRead).  This also makes it easier to provide
full specializations of the implementation.  "Full specializations"
mean that all the template parameters are fixed, so that the compiler
can actually compile the code.  This technique keeps your library's or
application's build times down, since kernels are already precompiled
for certain template parameter combinations.  It also improves
performance, since compilers have an easier time optimizing code in
shorter .cpp files.

## Building Kokkoskernels

### CMake
Following Kokkos style, all CMake options are of the form
````
KokkosKernels_ENABLE_OPTION
````
with options capitalized at the end. Almost all Kokkos Kernels options determine
whether ETI is used with a particular datatype, e.g.
````
-DKokkosKernels_INST_DOUBLE=On 
````
which does explicit instantation of all kernels for double type.
Kokkos Kernels derives most of its CXXFLAGS, C++ standard, architecture flags,
and other options from an installed (or in-tree) Kokkos package.
Tuning for a particular device or architecture is generally done through *Kokkos*
while tuning which kernels get instantiated is done through *Kokkos Kernels*.

Kokkos Kernels does supply flags for *asserting* properies of the linked Kokkos,
for example:
````
-DKokkosKernels_REQUIRE_DEVICES=CUDA
-DKokkosKernels_REQUIRE_OPTIONS=cuda_relocatable_device_code
````
This does *NOT* enable CUDA directly, but rather verifies that the underlying Kokkos supports
the desired option. If the underlying Kokkos was not built properly, CMake will crash
and tell you to re-build Kokkos. The values (unlike the option names) are not case-sensitive.
More details can be found in the [build instructions](BUILD.md) or [developer instructions](DEVELOPER.md).


### Spack

An alternative to manually building with the CMake is to use the Spack package manager.
To do so, download the `kokkos-spack` git repo and add to the package list:
````
spack repo add $path-to-kokkos-spack
````
A basic installation would be done as:
````
spack install kokkos-kernels
````
Spack allows options and and compilers to be tuned in the install command.
````
spack install kokkos-kernels@3.0 +double %gcc@7.3.0 +openmp
````
This example illustrates the three most common parameters to Spack:
* Variants: specified with, e.g. `+openmp`, this activates (or deactivates with, e.g. `~openmp`) certain options.
* Version:  immediately following `kokkos-kernels` the `@version` can specify a particular Kokkos Kernels to build
* Compiler: a default compiler will be chosen if not specified, but an exact compiler version can be given with the `%`option.

For a complete list of Kokkos Kernels options, run:
````
spack info kokkos-kernels
````

#### Tuning Kokkos Options
As discussed above in the CMake section, Kokkos Kernels inherits much of its configuration from the installed Kokkos.
Spack gives a mechanism for directly specifying Kokkos dependency options:
````
spack install kokkos-kernels ^kokkos@3.0+cuda+cuda_uvm
````
The carat `^` sepcifies an exact dependency configuration, which in this case activates CUDA and CUDA_UVM.
For a complete list of tunable Kokkos options, run
````
spack info kokkos
````

#### Setting up a development environment with Spack
Spack is generally most useful for installng packages to use.
If you want to install all *dependencies* of Kokkos Kernels first so that you can actively develop a given Kokkos Kernels source this can still be done. Go to the Kokkos Kernels source code folder and run:
````
spack diy -u cmake kokkos-kernels@{version} ...
````
specifying the exact version you want to develop and giving any spec options in `...`.
This creates a folder `spack-build` where you can `make`.

### Trilinos
For Trilinos builds with the Cuda backend and complex double enabled with ETI,
the cmake option below may need to be set to avoid Error 127 errors:
`CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS:BOOL=ON`

If the option above is not set, a warning will be issued during configuration:

"The CMake option CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS is either
undefined or OFF.  Please set
CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS:BOOL=ON when building with CUDA and
complex double enabled."


## Using Kokkoskernels Test Drivers 

In `perf_test` there are test drivers.

* `KokkosGraph_triangle.exe` : Triangle counting driver. 
* `KokkosSparse_spgemm.exe` : Sparse Matrix Sparse Matrix Multiply: 
* *****NOTE: KKMEM is outdated. Use default algorithm: KKSPGEMM = KKDEFAULT = DEFAULT****
* Or within the code:
*     kh.create_spgemm_handle(KokkosSparse::SPGEMM_KK);
* `KokkosSparse_spmv.exe` : Sparse matvec.
* `KokkosSparse_pcg.exe`  : CG method with Gauss Seidel as preconditioner.
* `KokkosGraph_color.exe` : Distance-1 Graph coloring 
* `KokkosKernels_MatrixConverter.exe` : given a matrix market format, converts it ".bin"
   binary format for fast input output readings, which can be read by other test drivers.
   
Please report bugs or performance issues to: https://github.com/kokkos/kokkos-kernels/issues

##### [LICENSE](https://github.com/kokkos/kokkos-kernels/blob/devel/LICENSE)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

KokkosKernels is licensed under standard 3-clause BSD terms of use.  For
specifics, please refer to the LICENSE file contained in the
repository or distribution.  Under the terms of Contract DE-NA0003525 with NTESS,
the U.S. Government retains certain rights in this software.
