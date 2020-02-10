#
# Disables across multiple builds on 'ride'
#

IF (ATDM_CMAKE_BUILD_TYPE STREQUAL "DEBUG")

  # Disable some expensive KokkosKernels tests in pure debug builds (#6464)
  ATDM_SET_ENABLE(KokkosKernels_sparse_serial_MPI_1_DISABLE ON)

ENDIF()

IF (Trilinos_ENABLE_DEBUG)

  # Disable Tempus tests that started timing out in debug builds when
  # Trilinos_ENABLE_DEBUG=ON was set PR #5970 (#6009)
  ATDM_SET_ENABLE(Tempus_BackwardEuler_MPI_1_DISABLE ON)
  ATDM_SET_ENABLE(Tempus_DIRK_ASA_MPI_1_DISABLE ON)
  ATDM_SET_ENABLE(Tempus_ExplicitRK_ASA_MPI_1_DISABLE ON)
  ATDM_SET_ENABLE(Tempus_HHTAlpha_MPI_1_DISABLE ON)
  ATDM_SET_ENABLE(Tempus_IMEX_RK_Combined_FSA_MPI_1_DISABLE ON)
  ATDM_SET_ENABLE(Tempus_IMEX_RK_Partitioned_Staggered_FSA_Partitioned_IMEX_RK_1st_Order_MPI_1_DISABLE ON)
  ATDM_SET_ENABLE(Tempus_IMEX_RK_Staggered_FSA_MPI_1_DISABLE ON)
  ATDM_SET_ENABLE(Tempus_Newmark_MPI_1_DISABLE ON)

ENDIF()

IF (ATDM_NODE_TYPE STREQUAL "CUDA")

  # Disable test that runs over 30 min currently (#2446)
  ATDM_SET_ENABLE(PanzerAdaptersSTK_MixedPoissonExample-ConvTest-Hex-Order-3_DISABLE ON)

  # Disable test that consistantly times out at 10 minutes (#3579)
  ATDM_SET_ENABLE(PanzerAdaptersSTK_CurlLaplacianExample-ConvTest-Quad-Order-4_DISABLE ON)

  # Disable randomly failing MueLu tests (#2311)
  ATDM_SET_ENABLE(MueLu_ParameterListInterpreterTpetra_MPI_1_DISABLE ON)
  ATDM_SET_ENABLE(MueLu_ParameterListInterpreterTpetraHeavy_MPI_1_DISABLE ON)

  IF (ATDM_CUDA_RDC)

    # Disable SEACAS tests (see #5784)
    ATDM_SET_ENABLE(pamgen_exodus_io_info_DISABLE ON)

  ENDIF()

ENDIF()