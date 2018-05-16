# #2410: This passes in the CUDA-RELESE-CUDA build so it is okay to disable in
# this DEBUG build.  This is likey related to the mix-lanauage compiler defect
# reported in #1208.
ATDM_SET_ENABLE(TeuchosNumerics_LAPACK_test_MPI_1_DISABLE ON)

# This test segfaults in the 'debug' builds on this system (#2466)
ATDM_SET_ENABLE(Belos_Tpetra_PseudoBlockCG_hb_test_MPI_4_DISABLE ON)

INCLUDE("${CMAKE_CURRENT_LIST_DIR}/CUDA_COMMON_TWEAKS.cmake")
