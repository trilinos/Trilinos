IF(CMAKE_PROJECT_NAME STREQUAL "Trilinos")
  set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Whether to build Serial backend")

  if(NOT ${Trilinos_ENABLE_OpenMP} STREQUAL "")
    set(Kokkos_ENABLE_OPENMP ${Trilinos_ENABLE_OpenMP} CACHE BOOL "Whether to build OpenMP backend")
  else()
    set(Kokkos_ENABLE_OPENMP                       OFF CACHE BOOL "Whether to build OpenMP backend")
  endif()

  if(NOT ${Kokkos_ENABLE_Pthread} STREQUAL "")
    set(Kokkos_ENABLE_THREADS ${Kokkos_ENABLE_Pthread} CACHE BOOL "Whether to build C++ threads backend")
  else()
    set(Kokkos_ENABLE_THREADS                   OFF CACHE BOOL "Whether to build C++ threads backend")
  endif()

  #message(FATAL_ERROR "Kokkos_ENABLE_CUDA: ${Kokkos_ENABLE_CUDA}")

  if(NOT ${Kokkos_ENABLE_CUDA} STREQUAL "")
    set(Kokkos_ENABLE_CUDA ${Kokkos_ENABLE_CUDA} CACHE BOOL "Whether to build CUDA backend")
  else()
    set(Kokkos_ENABLE_CUDA OFF)
    set(Kokkos_ENABLE_CUDA                OFF CACHE BOOL "Whether to build CUDA backend")
  endif()

  if(NOT ${Kokkos_ENABLE_HPX} STREQUAL "")
    set(Kokkos_ENABLE_HPX ${Kokkos_ENABLE_HPX} CACHE BOOL "Whether to build HPX backend")
  else()
    set(Kokkos_ENABLE_HPX OFF)
    set(Kokkos_ENABLE_HPX               OFF CACHE BOOL "Whether to build HPX backend")
  endif()

  if(NOT ${Kokkos_ENABLE_quadmath} STREQUAL "")
    set(Kokkos_ENABLE_LIBQUADMATH ${Kokkos_ENABLE_quadmath} CACHE BOOL "Whether to enable the LIBQUADMATH library")
  else()
    set(Kokkos_ENABLE_LIBQUADMATH                    OFF CACHE BOOL "Whether to enable the LIBQUADMATH library")
  endif()

  if(NOT ${Kokkos_ENABLE_DLlib} STREQUAL "")
    set(Kokkos_ENABLE_LIBDL ${Kokkos_ENABLE_DLlib} CACHE BOOL "Whether to enable the LIBDL library")
  else()
    set(Kokkos_ENABLE_LIBDL                 OFF CACHE BOOL "Whether to enable the LIBDL library")
  endif()

  set(Kokkos_ENABLE_COMPLEX_ALIGN OFF CACHE BOOL "Whether to align Kokkos::complex to 2*alignof(RealType)")

  set(Kokkos_ENABLE_TESTS OFF CACHE BOOL "Whether to build the unit tests")

  # FIXME_TRIBITS This should not be necessary if we define the correct targets
  tribits_package_postprocess()
ENDIF()
