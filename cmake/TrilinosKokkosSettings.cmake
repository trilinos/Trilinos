set(Kokkos_ENABLE_COMPLEX_ALIGN OFF CACHE BOOL "Whether to align Kokkos::complex to 2*alignof(RealType)")

# Enable serial backed for all builds
set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Whether to build Serial backend" FORCE)

# Trilinos_ENABLE_OpenMP -> Kokkos_ENABLE_OPENMP
if(NOT "${Trilinos_ENABLE_OpenMP}" STREQUAL "")
  set(Kokkos_ENABLE_OPENMP ${Trilinos_ENABLE_OpenMP} CACHE BOOL "Whether to build OpenMP backend" FORCE)
else()
  set(Kokkos_ENABLE_OPENMP OFF CACHE BOOL "Whether to build OpenMP backend")
endif()

# TPL_ENABLE_CUDA -> Kokkos_ENABLE_CUDA
if(NOT "${TPL_ENABLE_CUDA}" STREQUAL "")
  set(Kokkos_ENABLE_CUDA ${TPL_ENABLE_CUDA} CACHE BOOL "Whether to build CUDA backend" FORCE)
else()
  set(Kokkos_ENABLE_CUDA OFF CACHE BOOL "Whether to build CUDA backend")
endif()

# TPL_ENABLE_HPX -> Kokkos_ENABLE_HPX
if(NOT "${TPL_ENABLE_HPX}" STREQUAL "")
  set(Kokkos_ENABLE_HPX ${TPL_ENABLE_HPX} CACHE BOOL "Whether to build HPX backend" FORCE)
else()
  set(Kokkos_ENABLE_HPX OFF CACHE BOOL "Whether to build HPX backend")
endif()


# If only the Serial backend is enabled, turn off atomic bypass, unless the user set that option.
# These variables are undefined if the user did not set them.
if(Kokkos_ENABLE_SERIAL
  AND NOT Kokkos_ENABLE_OPENMP
  AND NOT Kokkos_ENABLE_THREADS
  AND NOT Kokkos_ENABLE_HPX
  AND NOT Kokkos_ENABLE_CUDA
  AND NOT Kokkos_ENABLE_HIP
  AND NOT Kokkos_ENABLE_SYCL
  AND NOT Kokkos_ENABLE_OPENACC
)
  if(NOT DEFINED Kokkos_ENABLE_ATOMICS_BYPASS)
    set(Kokkos_ENABLE_ATOMICS_BYPASS ON
      CACHE BOOL "Whether to make atomics non-atomic for non-threaded MPI-only use cases" FORCE
    )
  endif()
endif()
