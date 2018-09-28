#disable from issue #2410
ATDM_SET_ENABLE(TeuchosNumerics_LAPACK_test_MPI_1_DISABLE ON)

#disable from issue #2466
ATDM_SET_ENABLE(Belos_Tpetra_PseudoBlockCG_hb_test_MPI_4_DISABLE ON)

#disable failing test as in #3173
ATDM_SET_ENABLE(KokkosKernels_sparse_openmp_MPI_1_DISABLE ON)

INCLUDE("${CMAKE_CURRENT_LIST_DIR}/ALL_COMMON_TWEAKS.cmake")
