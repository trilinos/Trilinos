# #2410: STEQR() test fails on IBM Power systems with current TPL setup
ATDM_SET_ENABLE(TeuchosNumerics_DISABLE_STEQR_TEST ON)

# Disable timing out unit test in this debug build (#3336)
ATDM_SET_ENABLE(Kokkos_ContainersUnitTest_Serial_MPI_1_DISABLE ON)

# Disable some unit tests that run too slow in this DEBUG build (#2827)
ATDM_SET_CACHE(KokkosKernels_sparse_serial_MPI_1_EXTRA_ARGS
  "--gtest_filter=-serial.sparse_block_gauss_seidel_double_int_int_TestExecSpace:serial.sparse_block_gauss_seidel_double_int_size_t_TestExecSpace:serial.sparse_trsv_mv_double_int_int_LayoutLeft_TestExecSpace"
  CACHE STRING )
