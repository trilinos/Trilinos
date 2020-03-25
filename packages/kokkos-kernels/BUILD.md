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
That's it! There is no checking Kokkos preprocessor, compiler, or linker flags.
Kokkos propagates all the necesssary flags to your project.
This means not only is linking to Kokkos Kernels easy, but Kokkos Kernels itself can actually configure compiler and linker flags for *your* 
project. If building in-tree, there is no `find_package` and you link with `target_link_libraries(kokkoskernels)`.
In fact, you only ever need to link to Kokkos Kernels and not Kokkos!
Linking to Kokkos Kernels transitively provides Kokkos.


## Configuring CMake
A very basic installation is done with:
````
cmake ${srcdir} \
 -DCMAKE_CXX_COMPILER=g++ \
 -DCMAKE_INSTALL_PREFIX=${my_install_folder}
````
which builds and installs a default Kokkos Kernels when you run `make install`.
There are a few ETI (explicit template instantiation) and library options:
````
cmake ${srcdir} \
 -DCMAKE_CXX_COMPILER=g++ \
 -DCMAKE_INSTALL_PREFIX=${my_install_folder} \
 -DKokkosKernels_ENABLE_TPL_BLAS=ON
```` 
which activates the BLAS dependency.

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

# Kokkos Keyword Listing

## Device Backends
Device backends can be enabled by specifiying `-DKokkos_ENABLE_X`.

* KokkosKernels_ENABLE_EXAMPLES
  * Whether to build the example
  * BOOL Default: OFF
* KokkosKernels_ENABLE_EXPERIMENTAL
  * Whether to build experimental features
  * BOOL Default: OFF
* KokkosKernels_ENABLE_TESTS
  * Whether to build the tests
  * BOOL Default: OFF
* KokkosKernels_ENABLE_TPL_BLAS
  * Whether to link to a system BLAS
  * BOOL Default: OFF
* KokkosKernels_ENABLE_TPL_MKL
  * Whether to link to a system MKL
  * BOOL Default: OFF
* KokkosKernels_ENABLE_TPL_MAGMA
  * Whether to link to a system MAGMA
  * BOOL Default: OFF
* KokkosKernels_ENABLE_TPL_CUBLAS
  * Whether to link to the CUDA cuBLAS library
  * BOOL Default: ON if CUDA build, OFF otherwise
* KokkosKernels_ENABLE_TPL_CUSPARSE
  * Whether to link to the CUDA cuSPARSE library
  * BOOL Default: ON if CUDA build, OFF otherwise
* KokkosKernels_ETI_ONLY
  * Whether to restrict available kernels only to those that were explicitly pre-compiled with ETI
  * BOOL Default: OFF
* KokkosKernels_TEST_ETI_ONLY
  * Whether to restrict test kernels only to those that were explicitly pre-compiled with ETI
  * BOOL Default: ON

##### [LICENSE](https://github.com/kokkos/kokkos/blob/devel/LICENSE)

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Under the terms of Contract DE-NA0003525 with NTESS,
the U.S. Government retains certain rights in this software.
