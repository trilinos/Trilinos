![Kokkos](https://avatars2.githubusercontent.com/u/10199860?s=200&v=4)

# Installing and Using Kokkos Kernels

## Kokkos Kernels Philosophy
Kokkos Kernels provides a modern CMake style build system.
As C++ continues to develop for C++20 and beyond, CMake is likely to provide the most robust support
for C++.  Applications heavily leveraging Kokkos are strongly encouraged to use a CMake build system.

You can either use Kokkos Kernels as an installed package (encouraged) or use Kokkos in-tree in your project.
Modern CMake is exceedingly simple at a high-level (with the devil in the details).
Once Kokkos Kernels is installed In your `CMakeLists.txt` simply use:
````
find_package(KokkosKernels REQUIRED)
````
Then for every executable or library in your project:
````
target_link_libraries(myTarget Kokkos::kokkoskernels)
````
There is no checking Kokkos preprocessor, compiler, or linker flags.
CMake propagates all the necesssary flags to your project.
This means not only is linking to Kokkos Kernels easy, but Kokkos Kernels itself can actually configure compiler and linker flags for *your* 
project. If building in-tree, there is no `find_package`.
In fact, you only ever need to link to Kokkos Kernels and not Kokkos!
Linking to Kokkos Kernels transitively provides Kokkos.


## Configuring CMake
A very basic installation is done with:
````
cmake ${srcdir} \
 -DCMAKE_CXX_COMPILER=g++ \
 -DCMAKE_INSTALL_PREFIX=${my_install_folder}
 -DKokkos_ROOT=${kokkos_install_prefix}
````
which builds and installs a default Kokkos Kernels when you run `make install`.

There are also option Third-party Libraries (TPLs)
````
cmake ${srcdir} \
 -DCMAKE_CXX_COMPILER=g++ \
 -DCMAKE_INSTALL_PREFIX=${my_install_folder} \
 -DKokkosKernels_ENABLE_TPL_BLAS=ON
```` 
which, e.g. activates the BLAS dependency.
The full keyword listing is below.

## Spack
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
spack install kokkos-kernels@3.0 %gcc@7.3.0 +openmp
````
This example illustrates the three most common parameters to Spack:
* Variants: specified with, e.g. `+openmp`, this activates (or deactivates with, e.g. `~openmp`) certain options.
* Version:  immediately following `kokkos` the `@version` can specify a particular Kokkos to build
* Compiler: a default compiler will be chosen if not specified, but an exact compiler version can be given with the `%`option.

For a complete list of Kokkos Kernels options, run:
````
spack info kokkos-kernels
````

#### Spack Development
Spack currently installs packages to a location determined by a unique hash. This hash name is not really "human readable".
Generally, Spack usage should never really require you to reference the computer-generated unique install folder. 
If you must know, you can locate Spack Kokkos installations with:
````
spack find -p kokkos-kernels ...
````
where `...` is the unique spec identifying the particular Kokkos configuration and version.

A better way to use Spack for doing Kokkos development is the DIY feature of Spack.
If you wish to develop Kokkos Kernels itself, go to the Kokkos source folder:
````
spack diy -u cmake kokkos-kernels +diy ... 
````
where `...` is a Spack spec identifying the exact Kokkos Kernels configuration.
This then creates a `spack-build` directory where you can run `make`.

If you want more control on the underlying Kokkos, you can do:
````
spack diy -u cmake ${myproject}@${myversion} ... ^kokkos...
````
where the `...` are the specs for your project and the desired underlying Kokkos configuration.
Again, a `spack-build` directory will be created where you can run `make`.

Spack has a few idiosyncracies that make building outside of Spack annoying related to Spack forcing use of a compiler wrapper. This can be worked around by having a `-DSpack_WORKAROUND=On` given in your CMake. Then add the block of code to your CMakeLists.txt:

````
if (Spack_WORKAROUND)
 set(SPACK_CXX $ENV{SPACK_CXX})
 if(SPACK_CXX)
   set(CMAKE_CXX_COMPILER ${SPACK_CXX} CACHE STRING "the C++ compiler" FORCE)
   set(ENV{CXX} ${SPACK_CXX})
 endif()
endif()
````

# Kokkos Kernels CMake Option Listing

* BLAS_LIBRARIES: STRING
  * Optional override for the libraries that comprise TPL BLAS.
  * Default: None. Default common library names will be searched
* BLAS_LIBRARY_DIRS: STRING
  * Optional override for the library directories that comprise TPL BLAS.
  * Default: None. Default common library locations will be searched
* CUBLAS_LIBRARIES: STRING
  * Optional override for the libraries that comprise TPL CUBLAS.
  * Default: None. Default common library names will be searched
* CUBLAS_LIBRARY_DIRS: STRING
  * Optional override for the library directories that comprise TPL CUBLAS.
  * Default: None. Default common library locations will be searched
* CUSPARSE_LIBRARIES: STRING
  * Optional override for the libraries that comprise TPL CUSPARSE.
  * Default: None. Default common library names will be searched
* CUSPARSE_LIBRARY_DIRS: STRING
  * Optional override for the library directories that comprise TPL CUSPARSE.
  * Default: None. Default common library locations will be searched
* ARMPL_LIBRARIES: STRING
  * Optional override for the libraries that comprise TPL ARMPL.
  * Default: None. Default common library names will be searched
* ARMPL_LIBRARY_DIRS: STRING
  * Optional override for the library directories that comprise TPL ARMPL.
  * Default: None. Default common library locations will be searched
* KokkosKernels_BLAS_ROOT: PATH
  * Location of BLAS install root.
  * Default: None or the value of the environment variable BLAS_ROOT if set
* KokkosKernels_CUBLAS_ROOT: PATH
  * Location of CUBLAS install root.
  * Default: None or the value of the environment variable CUBLAS_ROOT if set
* KokkosKernels_CUSPARSE_ROOT: PATH
  * Location of CUSPARSE install root.
  * Default: None or the value of the environment variable CUSPARSE_ROOT if set
* KokkosKernels_ENABLE_EXAMPLES: BOOL
  * Whether to build examples.
  * Default: OFF
* KokkosKernels_ENABLE_EXPERIMENTAL: BOOL
  * Enable building and installation of experimental KokkosKernels features.
  * Default: OFF
* KokkosKernels_ENABLE_TESTS: BOOL
  * Whether to build tests.
  * Default: OFF
* KokkosKernels_ENABLE_PERFTESTS: BOOL
  * Whether to build performance tests.
  * Default: ON
* KokkosKernels_ENABLE_TESTS_AND_PERFSUITE: BOOL
  * Whether to build performance tests and suite.
  * Default: OFF
* KokkosKernels_ENABLE_DOCS: BOOL
  * Whether to build docs.
  * Default: OFF
* KokkosKernels_ENABLE_TPL_BLAS: BOOL
  * Whether to enable BLAS
  * Default: OFF
* KokkosKernels_ENABLE_TPL_CUBLAS: BOOL
  * Whether to enable CUBLAS
  * Default: ON if CUDA-enabled Kokkos, otherwise OFF
* KokkosKernels_ENABLE_TPL_CUSPARSE: BOOL
  * Whether to enable CUSPARSE
  * Default: ON if CUDA-enabled Kokkos, otherwise OFF
* KokkosKernels_ENABLE_TPL_LAPACK: BOOL
  * Whether to enable LAPACK
  * Default: ON if BLAS is enabled, otherwise OFF
* KokkosKernels_ENABLE_TPL_MAGMA: BOOL
  * Whether to enable MAGMA
  * Default: OFF
* KokkosKernels_ENABLE_TPL_MKL: BOOL
  * Whether to enable MKL
  * Default: OFF
* KokkosKernels_ENABLE_TPL_ARMPL: BOOL
  * Whether to enable ARMPL
  * Default: OFF
* KokkosKernels_ETI_ONLY: BOOL
  * Whether to restrict availability of kernels to ETI types only. Turning this on guarantees that kernels are never built inside of object files which simply call KokkosKernels functions.
  * Default: OFF
* KokkosKernels_INST_COMPLEX_DOUBLE: BOOL
  * Whether to pre instantiate kernels for the scalar type complex<double>.  Disabling this may increase build times.
  * Default: OFF or unless enabled during a Trilinos build with Trilinos_ENABLE_COMPLEX_DOUBLE.
* KokkosKernels_INST_COMPLEX_FLOAT: BOOL
  * Whether to pre instantiate kernels for the scalar type complex<float>.  Disabling this may increase build times.
  * Default: OFF or unless enabled during a Trilinos build with Trilinos_ENABLE_COMPLEX_FLOAT.
* KokkosKernels_INST_DOUBLE: BOOL
  * Whether to pre instantiate kernels for the scalar type double.  This option is KokkosKernels_INST_DOUBLE=ON by default.  Disabling this may increase build times.
  * Default: ON
* KokkosKernels_INST_EXECSPACE_OPENMP: BOOL
  * Whether to pre instantiate kernels for the execution space Kokkos::OpenMP.  Disabling this when Kokkos_ENABLE_OPENMP is enabled may increase build times.
  * Default: ON if Kokkos is OpenMP-enabled, OFF otherwise.
* KokkosKernels_INST_EXECSPACE_SERIAL: BOOL
  * Whether to build kernels for the execution space Kokkos::Serial.  If explicit template instantiation (ETI) is enabled in Trilinos, disabling this when Kokkos_ENABLE_SERIAL is enabled may increase build times.
  * Default: ON when Kokkos is Serial-enabled, OFF otherwise.
* KokkosKernels_INST_EXECSPACE_THREADS: BOOL
  * Whether to build kernels for the execution space Kokkos::Threads.  If explicit template instantiation (ETI) is enabled in Trilinos, disabling this when Kokkos_ENABLE_PTHREAD is enabled may increase build times.
  * Default: ON if Kokkos is Threads-enabled, OFF otherwise.
* KokkosKernels_INST_FLOAT: BOOL
  * Whether to pre instantiate kernels for the scalar type float.  Disabling this may increase build times.
  * Default: OFF or unless enabled during a Trilinos build with Trilinos_ENABLE_FLOAT.
* KokkosKernels_INST_LAYOUTLEFT: BOOL
  * Whether to pre instantiate kernels for the view layout LayoutLeft.  This option is KokkosKernels_INST_LAYOUTLEFT=ON by default.  Disabling this may increase build times.
  * Default: ON
* KokkosKernels_INST_LAYOUTRIGHT: BOOL
  * Whether to pre instantiate kernels for the view layout LayoutRight.  This option is KokkosKernels_INST_LAYOUTRIGHT=OFF by default.  Disabling this may increase build times.
  * Default: OFF
* KokkosKernels_INST_MEMSPACE_HOSTSPACE: BOOL
  * Whether to pre instantiate kernels for the memory space Kokkos::HostSpace.  Disabling this when one of the Host execution spaces is enabled may increase build times.
  * Default: ON
* KokkosKernels_INST_OFFSET_INT: BOOL
  * Whether to pre instantiate kernels for the offset type int. This option is KokkosKernels_INST_OFFSET_INT=ON by default.
  * Default: ON
* KokkosKernels_INST_OFFSET_SIZE_T: BOOL
  * Whether to pre instantiate kernels for the offset type size_t.  This option is KokkosKernels_INST_OFFSET_SIZE_T=OFF by default.
  * Default: ON
* KokkosKernels_INST_ORDINAL_INT: BOOL
  * Whether to pre instantiate kernels for the ordinal type int.  This option is KokkosKernels_INST_ORDINAL_INT=ON by default.
  * Default: ON
* KokkosKernels_INST_ORDINAL_INT64_T: BOOL
  * Whether to pre instantiate kernels for the ordinal type int64_t.  This option is KokkosKernels_INST_ORDINAL_INT64_T=OFF by default.
  * Default: OFF
* KokkosKernels_LAPACK_ROOT: PATH
  * Location of LAPACK install root.
  * Default: None or the value of the environment variable LAPACK_ROOT if set
* KokkosKernels_LINALG_OPT_LEVEL: BOOL **DEPRECATED**
  * Optimization level for KokkosKernels computational kernels: a nonnegative integer.  Higher levels result in better performance that is more uniform for corner cases, but increase build time and library size.  The default value is 1, which should give performance within ten percent of optimal on most platforms, for most problems.
  * Default: 1
* KokkosKernels_MAGMA_ROOT: PATH
  * Location of MAGMA install root.
  * Default: None or the value of the environment variable MAGMA_ROOT if set
* KokkosKernels_MKL_ROOT: PATH
  * Location of MKL install root.
  * Default: None or the value of the environment variable MKL_ROOT if set
* KokkosKernels_NO_DEFAULT_CUDA_TPLS: BOOL
  * Whether CUDA TPLs should be enabled by default.
  * Default: OFF
* KokkosKernels_TEST_ETI_ONLY: BOOL
  * Whether to restrict testing to ETI types.
  * Default: ON
* LAPACK_LIBRARIES: STRING
  * Optional override for the libraries that comprise TPL LAPACK.
  * Default: None. Default common library names will be searched
* LAPACK_LIBRARY_DIRS: STRING
  * Optional override for the library directories that comprise TPL LAPACK.
  * Default: None. Default common library locations will be searched
* MKL_LIBRARIES: STRING
  * Optional override for the libraries that comprise TPL MKL.
  * Default: None. Default common library names will be searched
* MKL_LIBRARY_DIRS: STRING
  * Optional override for the library directories that comprise TPL MKL.
  * Default: None. Default common library locations will be searched

##### [LICENSE](https://github.com/kokkos/kokkos/blob/devel/LICENSE)

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Under the terms of Contract DE-NA0003525 with NTESS,
the U.S. Government retains certain rights in this software.
