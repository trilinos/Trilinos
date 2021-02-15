#
# Set up to limit running on GPUs
#

ATDM_SET_CACHE(Trilinos_AUTOGENERATE_TEST_RESOURCE_FILE ON CACHE BOOL)
ATDM_SET_CACHE(Trilinos_CUDA_NUM_GPUS 2 CACHE STRING)
ATDM_SET_CACHE(Trilinos_CUDA_SLOTS_PER_GPU 2 CACHE STRING)

#
# Disables across multiple builds on 'weaver'
#
ATDM_SET_ENABLE(SEACASIoss_structured_cgns_assembly_copy_DISABLE ON)
ATDM_SET_ENABLE(SEACASIoss_structured_cgns_assembly_copy_fpp_DISABLE ON)
ATDM_SET_ENABLE(Intrepid2_unit-test_Discretization_Basis_HierarchicalBases_Hierarchical_Basis_Tests_MPI_1_DISABLE ON)

# Disable known failures for SPARC Trilinos configuration (#3632)
ATDM_SET_ENABLE(PanzerAdaptersIOSS_tIOSSConnManager2_MPI_2_DISABLE ON)
ATDM_SET_ENABLE(PanzerAdaptersIOSS_tIOSSConnManager3_MPI_3_DISABLE ON)

# Disable randomly timing out test in all 'weaver' builds (#6463)
ATDM_SET_ENABLE(Teko_testdriver_tpetra_MPI_4_DISABLE ON)

IF (ATDM_CMAKE_BUILD_TYPE STREQUAL "DEBUG")

  # Disable some expensive KokkosKernels tests in pure debug builds (#6464)
  ATDM_SET_ENABLE(KokkosKernels_sparse_serial_MPI_1_DISABLE ON)

ENDIF()

IF (Trilinos_ENABLE_DEBUG)

  # STEQR() test fails on IBM Power systems with current TPL setup (#2410, #6166)
  ATDM_SET_ENABLE(TeuchosNumerics_DISABLE_STEQR_TEST ON)

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

  # This test fails consistently (#2751)
  ATDM_SET_ENABLE(PanzerAdaptersSTK_MixedPoissonExample-ConvTest-Hex-Order-3_DISABLE ON)

  # Disable known falure for ROL CUDA builds (#3543)
  ATDM_SET_ENABLE(ROL_test_elementwise_TpetraMultiVector_MPI_4_DISABLE ON)

  # Disable known failure (#6329)
  ATDM_SET_ENABLE(ROL_NonlinearProblemTest_MPI_4_DISABLE ON)

  # Disable known randomly timing out test (#7090)
  ATDM_SET_ENABLE(MueLu_ParameterListInterpreterTpetra_MPI_1_DISABLE ON)

  IF (ATDM_CUDA_RDC)

    # Disable the build of SEACAS 'explore' for all cuda+rdc builds for now (#6008)
    ATDM_SET_ENABLE(explore_EXE_DISABLE ON)

  ENDIF()

ENDIF()
