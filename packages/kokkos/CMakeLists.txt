include(CMakeListsKokkos.txt)

if(NOT ${Trilinos_ENABLE_OpenMP} STREQUAL "")
  set(Kokkos_ENABLE_OPENMP ${Trilinos_ENABLE_OPENMP} CACHE BOOL "Whether to build OpenMP backend" FORCE)
endif()

if(NOT ${TPL_ENABLE_CUDA} STREQUAL "")
  set(Kokkos_ENABLE_CUDA ${TPL_ENABLE_CUDA} CACHE BOOL "Whether to build Cuda backend" FORCE)
endif()

if(NOT ${TPL_ENABLE_HPX} STREQUAL "")
  set(Kokkos_ENABLE_HPX ${TPL_ENABLE_HPX} CACHE BOOL "Whether to build HPX backend" FORCE)
endif()

if(NOT ${TPL_ENABLE_quadmath} STREQUAL "")
  set(Kokkos_ENABLE_LIQUADMATH ${TPL_ENABLE_quadmath} CACHE BOOL "Whether to enable the LIBQUADMATH library" FORCE)
endif()

if(NOT ${TPL_ENABLE_DLlib} STREQUAL "")
  set(Kokkos_ENABLE_LIBDL ${TPL_ENABLE_DLlib} CACHE BOOL "Whether to enable the LIBDL library" FORCE)
endif()

set(Kokkoe_ENABLE_COMPLEX_ALIGN OFF CACHE BOOL "Whether to align Kokkos::complex to 2*alignof(RealType)" FORCE)

add_library(Kokkos::all_libs ALIAS kokkos)
tribits_package_postprocess()
file(GLOB KOKKOS_CMAKE_FILES ${Kokkos_BINARY_DIR}/*.cmake)
file(COPY ${KOKKOS_CMAKE_FILES} DESTINATION ${Trilinos_BINARY_DIR}/cmake_packages/kokkos/)
