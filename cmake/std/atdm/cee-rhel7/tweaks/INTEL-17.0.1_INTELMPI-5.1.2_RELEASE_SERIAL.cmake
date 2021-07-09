# Disable test that does not link due to HDF5 errors only for Intel env (#3891)
ATDM_SET_ENABLE(SEACASIoss_Utst_structured_decomp_EXE_DISABLE ON)
ATDM_SET_ENABLE(SEACASIoss_Utst_structured_decomp_MPI_1_DISABLE ON)

INCLUDE("${CMAKE_CURRENT_LIST_DIR}/ALL_COMMON_TWEAKS.cmake")
