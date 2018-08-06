# #2410: This passes in the GNU-RELESE-OPENMP build so it is okay to disable
# in this DEBUG build.  This is likey related to the mix-lanauage compiler
# defect reported in #1208.
ATDM_SET_ENABLE(TeuchosNumerics_LAPACK_test_MPI_1_DISABLE ON)

# Disable test that runs over 30 min currently (#2446)
ATDM_SET_ENABLE(PanzerAdaptersSTK_MixedPoissonExample-ConvTest-Hex-Order-3_DISABLE ON)

# This test segfaults in the 'debug' builds on this system (#2466)
ATDM_SET_ENABLE(Belos_Tpetra_PseudoBlockCG_hb_test_MPI_4_DISABLE ON)

# Disable some unit tests that run too slow in this DEBUG build (#2827)
ATDM_SET_CACHE(KokkosContainers_UnitTest_Serial_MPI_1_EXTRA_ARGS
  "--gtest_filter=-serial.bitset:serial.scatterview"
  CACHE STRING )
ATDM_SET_CACHE(KokkosContainers_UnitTest_OpenMP_MPI_1_EXTRA_ARGS
  "--gtest_filter=-openmp.UnorderedMap_failed_insert"
  CACHE STRING )
ATDM_SET_CACHE(KokkosKernels_graph_serial_MPI_1_EXTRA_ARGS
  "--gtest_filter=-serial.graph_graph_color_d2_double_int_int_TestExecSpace"
  CACHE STRING )
ATDM_SET_CACHE(KokkosKernels_sparse_openmp_MPI_1_EXTRA_ARGS
  "--gtest_filter=-openmp.sparse_block_gauss_seidel_double_int_int_TestExecSpace:openmp.sparse_trsv_mv_double_int_int_LayoutLeft_TestExecSpace"
  CACHE STRING )
ATDM_SET_CACHE(KokkosKernels_sparse_serial_MPI_1_EXTRA_ARGS
  "--gtest_filter=-serial.sparse_block_gauss_seidel_double_int_int_TestExecSpace:serial.sparse_block_gauss_seidel_double_int_size_t_TestExecSpace:serial.sparse_trsv_mv_double_int_int_LayoutLeft_TestExecSpace"
  CACHE STRING )

# Disable entire test that is timing out (or nearly timing out) at 10 minutes
# in debug-openmp build (#3168)
ATDM_SET_ENABLE(KokkosKernels_sparse_openmp_MPI_1_DISABLE ON)

# Disable randomly failing tests for this build (#2920)
ATDM_SET_ENABLE(Belos_pseudo_stochastic_pcg_hb_0_MPI_4_DISABLE ON)
ATDM_SET_ENABLE(Belos_pseudo_stochastic_pcg_hb_1_MPI_4_DISABLE ON)

INCLUDE("${CMAKE_CURRENT_LIST_DIR}/ALL_COMMON_TWEAKS.cmake")
