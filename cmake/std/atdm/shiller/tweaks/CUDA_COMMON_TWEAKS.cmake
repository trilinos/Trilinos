# Disable test that runs over 30 min currently (#2446)
ATDM_SET_ENABLE(PanzerAdaptersSTK_MixedPoissonExample-ConvTest-Hex-Order-3_DISABLE ON)

# Disable test that times out after 10 minutes (#2455)
ATDM_SET_ENABLE(Anasazi_Epetra_BlockDavidson_auxtest_MPI_4_DISABLE ON)

# Disable test takes a long time to complete for some reason (#2455)
ATDM_SET_ENABLE(Anasazi_Epetra_LOBPCG_auxtest_MPI_4_DISABLE ON)

# Disable expensive test that started timing out (#2919)
ATDM_SET_ENABLE(Belos_rcg_hb_MPI_4_DISABLE ON)
