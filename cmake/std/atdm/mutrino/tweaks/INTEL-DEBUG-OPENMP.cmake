# Disable test that is randomly failing in this build (#2474)
ATDM_SET_ENABLE(Piro_MatrixFreeDecorator_UnitTests_MPI_4_DISABLE ON)

# Disable MueLu tests having file I/O problems on mutrinos (#2839)
ATDM_SET_ENABLE(MueLu_CreateOperatorTpetra_MPI_1_DISABLE ON)
ATDM_SET_ENABLE(MueLu_CreateOperatorTpetra_MPI_4_DISABLE ON)
ATDM_SET_ENABLE(MueLu_ParameterListInterpreterTpetra_MPI_1_DISABLE ON)
ATDM_SET_ENABLE(MueLu_ParameterListInterpreterTpetra_MPI_4_DISABLE ON)
ATDM_SET_ENABLE(MueLu_ParameterListInterpreterTpetraHeavy_MPI_1_DISABLE ON)
ATDM_SET_ENABLE(MueLu_ParameterListInterpreterTpetraHeavy_MPI_4_DISABLE ON)

# Disable SEACAS test that fails on mutrino (#2815)
ATDM_SET_ENABLE(SEACASExodus_exodus_unit_tests_DISABLE ON)
