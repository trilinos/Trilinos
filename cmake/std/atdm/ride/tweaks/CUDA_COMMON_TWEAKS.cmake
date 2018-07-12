# Disable test that runs over 30 min currently (#2446)
ATDM_SET_ENABLE(PanzerAdaptersSTK_MixedPoissonExample-ConvTest-Hex-Order-3_DISABLE ON)

# Disable randomly failing MueLu tests (#2311)
ATDM_SET_ENABLE(MueLu_ParameterListInterpreterTpetra_MPI_1_DISABLE ON)
ATDM_SET_ENABLE(MueLu_ParameterListInterpreterTpetraHeavy_MPI_1_DISABLE ON)
