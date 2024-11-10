tribits_package_define_dependencies(
  LIB_REQUIRED_PACKAGES
    Kokkos KokkosKernels
  LIB_OPTIONAL_TPLS
    MPI CUDA
  TEST_REQUIRED_PACKAGES
    Gtest
  )
