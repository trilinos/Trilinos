# Disable tests for all builds for this system

# Disable SEACAS tests that get messed up due to extra STDERR output on
# 'mutrino' that does not occur on other platforms (see #3183)
ATDM_SET_ENABLE(SEACASAprepro_aprepro_array_test_DISABLE ON)
ATDM_SET_ENABLE(SEACASAprepro_aprepro_command_line_include_test_DISABLE ON)
ATDM_SET_ENABLE(SEACASAprepro_aprepro_command_line_vars_test_DISABLE ON)
ATDM_SET_ENABLE(SEACASAprepro_aprepro_test_dump_reread_DISABLE ON)
ATDM_SET_ENABLE(SEACASAprepro_aprepro_unit_test_DISABLE ON)
ATDM_SET_ENABLE(SEACASAprepro_lib_aprepro_lib_array_test_DISABLE ON)
ATDM_SET_ENABLE(SEACASAprepro_lib_aprepro_lib_unit_test_DISABLE ON)
ATDM_SET_ENABLE(SEACASExodus_exodus_unit_tests_nc5_env_DISABLE ON)

# Disable seaas tests with strange Not Run with a not found for a command (see
# #3496)
ATDM_SET_ENABLE(SEACASAprepro_aprepro_test_exodus_DISABLE ON)
ATDM_SET_ENABLE(SEACASIoss_exodus32_to_exodus32_DISABLE ON)
ATDM_SET_ENABLE(SEACASIoss_exodus32_to_exodus32_pnetcdf_DISABLE ON)
ATDM_SET_ENABLE(SEACASIoss_exodus32_to_exodus64_DISABLE ON)

# Disable muelu tests that fail to build due to
# '"Kokkos::Compat" has no member "KokkosSerialWrapperNode"'
#ATDM_SET_ENABLE(MueLu_Maxwell3D-Tpetra_MPI_4_DISABLE ON)
#ATDM_SET_ENABLE(MueLu_Maxwell3D_EXE_DISABLE ON)

#message("ATDM_NODE_TYPE=${ATDM_NODE_TYPE}")
#message("ATDM_CONFIG_KOKKOS_ARCH=$ENV{ATDM_CONFIG_KOKKOS_ARCH}")
##message("Kokkos_ARCH_HSW=${Kokkos_ARCH_HSW}")

IF (ATDM_NODE_TYPE STREQUAL "OPENMP")
  # Disable tests for all OpenMP builds for this system
  IF ("$ENV{ATDM_CONFIG_KOKKOS_ARCH}" STREQUAL "HSW")
    # Disable SEACAS test that fails on mutrino (#2815)
    ATDM_SET_ENABLE(SEACASExodus_exodus_unit_tests_DISABLE ON)
  ENDIF()
ENDIF()

# Disable all ROL tests for ats1+hsw+intel-19+openmp (#7470)
IF (
    (ATDM_COMPILER MATCHES "INTEL-19") AND
    (ATDM_KOKKOS_ARCH STREQUAL "HSW") AND
    (ATDM_NODE_TYPE STREQUAL "OPENMP")
  )
  ATDM_SET_ENABLE(ROL_SKIP_CTEST_ADD_TEST ON)
ENDIF()
