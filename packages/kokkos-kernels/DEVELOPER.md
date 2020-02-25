![Kokkos](https://avatars2.githubusercontent.com/u/10199860?s=200&v=4)

# Developing Kokkos Kernels

This document contains build instructions for developers.
For build system details for users, refer to the [build instructions](BUILD.md).

## Build System

Kokkos Kernels uses CMake to configure, build, and install.
For a full overview of Modern CMake practices used by Kokkos projects,
consult the DEVELOPER.md documentation in the Kokkos repo.
Here we review a few extra details specific to Kokkos Kernels.

## Kokkos Kernels Options
Kokkos Kernels provide a set of user visible options `KokkosKernels_X` prefixed by a camel case name and followed by an all uppercase option name.
These camel case names are CMake cache variables, which can be given by `-DKokkosKernels_X=VALUE` on the command line.
These variables appear with docstrings in the `CMakeCache.txt` or can be edited with `ccmake`.
Because of the way CMake caches persist, best practice is to have a minimal set of cache variables and rely on regular variables to implement most of the logic.
Variables appearing in all-caps are CMake regular variables, not appearing in the cache.
To add a user-visible option, we use the provided function, e.g.
````
KOKKOSKERNELS_ADD_OPTION(
  ENABLE_THING
  ON
  BOOL
  "Enable a thing in the configuration"
)
````
When running CMake, this will look for a user-provided cache variable `KokkosKernels_ENABLE_THING`.
If not found, it sets to a default value `ON` - which is a `BOOL`. 
`STRING` is the other common option type. The function additionally creates a regular CMake variable `KOKKOSKERNELS_ENABLE_THING` with the given value.

## Transitive Dependence on Kokkos

Kokkos Kernels should never need to select its own `CXXFLAGS` or include directories.
The Kokkos library contains all the logic for selecting C++ standard, architecture, or other flags.
Through the magic of CMake targets and properties, the following line transitively adds these to Kokkos Kernels:
````
TARGET_LINK_LIBRARIES(kokkoskernels PUBLIC Kokkos::kokkos)
````

## Third-party Libraries (TPLs)

Kokkos Kernels can use BLAS (or MKL) or if using CUDA also cuBLAS or cuSPARSE.
An extra function is available for detecting whether the user wants to enable TPLs.
````
KOKKOSKERNELS_TPL_OPTION(
  BLAS
  OFF
  "Enable BLAS"
)
````
This creates a user-visible option `KokkosKernels_ENABLE_TPL_BLAS` with the default `OFF`.
Here the value is always a `BOOL`. The developer should add the necessary custom logic for that TPL, e.g.
````
IF (KOKKOSKERNELS_ENABLE_TPL_BLAS)
  FIND_PACKAGE(BLAS)
  ...
ENDIF()
````


