#
# Disables across multiple builds on 'ride'
#

IF (Trilinos_ENABLE_DEBUG)

  # Disable Tempus tests that started timing out in debug builds when
  # Trilinos_ENABLE_DEBUG=ON was set PR #5970 (#6009)
  ATDM_SET_ENABLE(Tempus_BackwardEuler_MPI_1_DISABLE ON)
  ATDM_SET_ENABLE(Tempus_DIRK_ASA_MPI_1_DISABLE ON)
  ATDM_SET_ENABLE(Tempus_ExplicitRK_ASA_MPI_1_DISABLE ON)
  ATDM_SET_ENABLE(Tempus_HHTAlpha_MPI_1_DISABLE ON)
  ATDM_SET_ENABLE(Tempus_IMEX_RK_Combined_FSA_MPI_1_DISABLE ON)
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

  #
  # Disable tests already failing in the non-RDC build
  #

  IF (ATDM_CUDA_RDC)

    # Disable ROL tests (see #3543)
    ATDM_SET_ENABLE(ROL_example_PDE-OPT_0ld_adv-diff-react_example_01_MPI_4_DISABLE ON)
    ATDM_SET_ENABLE(ROL_example_PDE-OPT_0ld_adv-diff-react_example_02_MPI_4_DISABLE ON)
    ATDM_SET_ENABLE(ROL_example_PDE-OPT_0ld_poisson_example_01_MPI_4_DISABLE ON)
    ATDM_SET_ENABLE(ROL_example_PDE-OPT_0ld_stefan-boltzmann_example_03_MPI_4_DISABLE ON)
    ATDM_SET_ENABLE(ROL_example_PDE-OPT_navier-stokes_example_01_MPI_4_DISABLE ON)
    ATDM_SET_ENABLE(ROL_example_PDE-OPT_navier-stokes_example_02_MPI_4_DISABLE ON)
    ATDM_SET_ENABLE(ROL_example_PDE-OPT_nonlinear-elliptic_example_01_MPI_4_DISABLE ON)
    ATDM_SET_ENABLE(ROL_example_PDE-OPT_nonlinear-elliptic_example_02_MPI_4_DISABLE ON)
    ATDM_SET_ENABLE(ROL_example_PDE-OPT_obstacle_example_01_MPI_4_DISABLE ON)
    ATDM_SET_ENABLE(ROL_example_PDE-OPT_stefan-boltzmann_example_01_MPI_4_DISABLE ON)
    ATDM_SET_ENABLE(ROL_example_PDE-OPT_stefan-boltzmann_example_03_MPI_4_DISABLE ON)
    ATDM_SET_ENABLE(ROL_example_PDE-OPT_topo-opt_poisson_example_01_MPI_4_DISABLE ON)
    ATDM_SET_ENABLE(ROL_test_elementwise_TpetraMultiVector_MPI_4_DISABLE ON)

    # Disable Zoltan tests (see #3749)
    ATDM_SET_ENABLE(TrilinosCouplings_Example_Maxwell_MueLu_MPI_1_DISABLE ON)
    ATDM_SET_ENABLE(TrilinosCouplings_Example_Maxwell_MueLu_MPI_4_DISABLE ON)

    # Disable Zoltan tests (see #4042)
    ATDM_SET_ENABLE(Zoltan_ch_ewgt_zoltan_parallel_DISABLE ON)
    ATDM_SET_ENABLE(Zoltan_ch_grid20x19_zoltan_parallel_DISABLE ON)
    ATDM_SET_ENABLE(Zoltan_ch_nograph_zoltan_parallel_DISABLE ON)
    ATDM_SET_ENABLE(Zoltan_ch_simple_zoltan_parallel_DISABLE ON)

    # Disable SEACAS tests (see #5784)
    ATDM_SET_ENABLE(pamgen_exodus_io_info_DISABLE ON)

  ENDIF()

ENDIF()