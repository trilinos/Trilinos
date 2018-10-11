# Disable test that is randomly failing in this build (#2474)
ATDM_SET_ENABLE(Piro_MatrixFreeDecorator_UnitTests_MPI_4_DISABLE ON)

# Disable MueLu tests having file I/O problems on mutrino (#2839)
ATDM_SET_ENABLE(MueLu_CreateOperatorTpetra_MPI_1_DISABLE ON)
ATDM_SET_ENABLE(MueLu_CreateOperatorTpetra_MPI_4_DISABLE ON)
ATDM_SET_ENABLE(MueLu_ParameterListInterpreterTpetra_MPI_1_DISABLE ON)
ATDM_SET_ENABLE(MueLu_ParameterListInterpreterTpetra_MPI_4_DISABLE ON)
ATDM_SET_ENABLE(MueLu_ParameterListInterpreterTpetraHeavy_MPI_1_DISABLE ON)
ATDM_SET_ENABLE(MueLu_ParameterListInterpreterTpetraHeavy_MPI_4_DISABLE ON)

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
