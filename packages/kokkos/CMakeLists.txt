include(CMakeListsKokkos.txt)

if(NOT ${Trilinos_ENABLE_OpenMP} STREQUAL "")
  kokkos_set_option(ENABLE_OPENMP ${Trilinos_ENABLE_OpenMP})
endif()

if(NOT ${TPL_ENABLE_CUDA} STREQUAL "")
  kokkos_set_option(ENABLE_CUDA ${TPL_ENABLE_CUDA})
endif()

kokkos_set_option(ENABLE_COMPLEX_ALIGN OFF)

add_library(Kokkos::all_libs ALIAS kokkos)
tribits_package_postprocess()
file(GLOB KOKKOS_CMAKE_FILES ${Kokkos_BINARY_DIR}/*.cmake)
file(COPY ${KOKKOS_CMAKE_FILES} DESTINATION ${Trilinos_BINARY_DIR}/cmake_packages/kokkos/)
