# Disable known failures for SPARC Trilinos configuration (#3632)
ATDM_SET_ENABLE(PanzerAdaptersIOSS_tIOSSConnManager2_MPI_2_DISABLE ON)
ATDM_SET_ENABLE(PanzerAdaptersIOSS_tIOSSConnManager3_MPI_3_DISABLE ON)

IF (ATDM_NODE_TYPE STREQUAL "CUDA")

  # This test fails consistently (#2751)
  ATDM_SET_ENABLE(PanzerAdaptersSTK_MixedPoissonExample-ConvTest-Hex-Order-3_DISABLE ON)

  # Disable known falure for ROL CUDA builds (#3543)
  ATDM_SET_ENABLE(ROL_test_elementwise_TpetraMultiVector_MPI_4_DISABLE ON)

ENDIF()
