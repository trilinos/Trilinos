#disable from issue #2410
ATDM_SET_ENABLE(TeuchosNumerics_LAPACK_test_MPI_1_DISABLE ON)

#disable from issue #2466
ATDM_SET_ENABLE(Belos_Tpetra_PseudoBlockCG_hb_test_MPI_4_DISABLE ON)

# Disable timing out unit test in this debug build (#3336)
ATDM_SET_ENABLE(KokkosContainers_UnitTest_Serial_MPI_1_DISABLE ON)

# Disable some unit tests that run too slow in this DEBUG build (#2827)
ATDM_SET_CACHE(KokkosKernels_sparse_serial_MPI_1_EXTRA_ARGS
  "--gtest_filter=-serial.sparse_block_gauss_seidel_double_int_int_TestExecSpace:serial.sparse_block_gauss_seidel_double_int_size_t_TestExecSpace:serial.sparse_trsv_mv_double_int_int_LayoutLeft_TestExecSpace"
  CACHE STRING )

INCLUDE("${CMAKE_CURRENT_LIST_DIR}/CUDA_COMMON_TWEAKS.cmake")
INCLUDE("${CMAKE_CURRENT_LIST_DIR}/ALL_COMMON_TWEAKS.cmake")
