![Kokkos](https://avatars2.githubusercontent.com/u/10199860?s=200&v=4)

# How to Setup CMake for Your Project

## Kokkos Kernels Philosophy
Kokkos Kernels provides a modern CMake style build system.
As C++ continues to develop for C++20 and beyond, CMake is likely to provide the most robust support
for C++.  Applications heavily leveraging Kokkos are strongly encouraged to use a CMake build system.

Contained in this folder are two implementations of the same project.
One uses `find_package` and links to installed Kokkos libraries.
The other builds Kokkos and Kokkos Kernels directly as part of the project,
which would occur, e.g., if they were submodules.

## Finding Installed Packages

CMake installs file called `<Package>Config.cmake` in the chosen install folder for a given `<PACKAGE>`.
For Kokkos, it installs `KokkosConfig.cmake` in `<prefix>/lib/cmake/Kokkos` and for Kokkos Kernels it installs a `KokkosKernelsConfig.cmake` in `<prefix>/lib/cmake/KokkosKernels`.

For `find_package` calls in the `CMakeLists.txt` to succeed, you must point CMake to the correct installation folder.
Prior to CMake 3.12, you would use:

````
cmake <my_project_source> -DKokkos_DIR=<kokkos_prefix>/lib/cmake/Kokkos
````
From 3.12 on, you can use:
````
cmake <my_project_source> -DKokkos_ROOT=<kokkos_prefix>
````
For projects depending on KokkosKernels, equivalent `KokkosKernels_ROOT` or `KokkosKernels_DIR` flags are needed.


