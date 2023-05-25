if(NOT ${Trilinos_ENABLE_OpenMP} STREQUAL "")
  set(Kokkos_ENABLE_OPENMP ${Trilinos_ENABLE_OpenMP} CACHE BOOL "Whether to build OpenMP backend" FORCE)
endif()

if(NOT ${TPL_ENABLE_CUDA} STREQUAL "")
  set(Kokkos_ENABLE_CUDA ${TPL_ENABLE_CUDA} CACHE BOOL "Whether to build CUDA backend" FORCE)
endif()

if(NOT ${TPL_ENABLE_HPX} STREQUAL "")
  set(Kokkos_ENABLE_HPX ${TPL_ENABLE_HPX} CACHE BOOL "Whether to build HPX backend" FORCE)
endif()

if(NOT ${TPL_ENABLE_quadmath} STREQUAL "")
  set(Kokkos_ENABLE_LIBQUADMATH ${TPL_ENABLE_OpenMP} CACHE BOOL "Whether to enable the LIBQUADMATH library" FORCE)
endif()

if(NOT ${TPL_ENABLE_DLlib} STREQUAL "")
  set(Kokkos_ENABLE_LIBDL ${TPL_ENABLE_DLlib} CACHE BOOL "Whether to enable the LIBDL library" FORCE)
endif()

set(Kokkos_ENABLE_COMPLEX_ALIGN OFF CACHE BOOL "Whether to build OpenMP backend" FORCE)

kokkos_set_option(ENABLE_COMPLEX_ALIGN OFF FORCE)
