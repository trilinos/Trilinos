#
# Disables across multiple builds on 'ats2'
#

IF (ATDM_CMAKE_BUILD_TYPE STREQUAL "DEBUG")

  # Disable some expensive KokkosKernels tests in pure debug builds (#6464)
  ATDM_SET_ENABLE(KokkosKernels_sparse_serial_MPI_1_DISABLE ON)

  # Disable test that just has problems with BinUtils interaction (#6821)
  ATDM_SET_ENABLE(TeuchosCore_show_stack_DISABLE ON)

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

  # Disble some timing out ROL tests (#6124)
  ATDM_SET_ENABLE(ROL_example_PDE-OPT_helmholtz_example_02_MPI_1_DISABLE ON)
  ATDM_SET_ENABLE(ROL_example_PDE-OPT_navier-stokes_example_01_MPI_4_DISABLE ON)

ENDIF()
