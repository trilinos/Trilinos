include(CMakeListsKokkos.txt)

kokkos_set_option(ENABLE_OPENMP ${Trilinos_ENABLE_OpenMP})
if(NOT ${Trilinos_ENABLE_OpenMP} STREQUAL "")
  kokkos_set_option(ENABLE_OPENMP ${Trilinos_ENABLE_OpenMP})
endif()

if(NOT ${TPL_ENABLE_CUDA} STREQUAL "")
  kokkos_set_option(ENABLE_CUDA ${TPL_ENABLE_CUDA})
endif()

if(NOT ${TPL_ENABLE_HPX} STREQUAL "")
  kokkos_set_option(ENABLE_HPX ${TPL_ENABLE_HPX})
endif()

if(NOT ${TPL_ENABLE_quadmath} STREQUAL "")
  kokkos_set_option(ENABLE_LIQUADMATH ${TPL_ENABLE_quadmath})
endif()

if(NOT ${TPL_ENABLE_DLlib} STREQUAL "")
  kokkos_set_option(ENABLE_LIBDL ${TPL_ENABLE_DLlib})
endif()

kokkos_set_option(ENABLE_COMPLEX_ALIGN OFF)

add_library(Kokkos::all_libs ALIAS kokkos)
tribits_package_postprocess()
file(GLOB KOKKOS_CMAKE_FILES ${Kokkos_BINARY_DIR}/*.cmake)
file(COPY ${KOKKOS_CMAKE_FILES} DESTINATION ${Trilinos_BINARY_DIR}/cmake_packages/kokkos/)
