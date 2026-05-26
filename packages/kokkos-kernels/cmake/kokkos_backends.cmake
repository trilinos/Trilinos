# Kokkos only defines the variables if the backends are ON
# Define aux variables that exist as on/off
macro(check_kokkos_backend BE)
  if(Kokkos_ENABLE_${BE})
    set(KOKKOS_ENABLE_${BE} ON)
  else()
    set(KOKKOS_ENABLE_${BE} OFF)
  endif()
  set(KOKKOSKERNELS_INST_EXECSPACE_${BE}_DEFAULT ${KOKKOS_ENABLE_${BE}})
endmacro()

check_kokkos_backend(SERIAL)
check_kokkos_backend(THREADS)
check_kokkos_backend(OPENMP)
check_kokkos_backend(OPENMPTARGET)
check_kokkos_backend(CUDA)
check_kokkos_backend(HIP)
check_kokkos_backend(SYCL)
