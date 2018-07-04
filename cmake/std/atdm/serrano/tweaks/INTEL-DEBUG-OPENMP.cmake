# Disable test that is failing or timing out in this build (see #2751)
ATDM_SET_ENABLE(PanzerAdaptersSTK_MixedPoissonExample-ConvTest-Hex-Order-3_DISABLE ON)

# Disable test that is randomly failing in this build (#2474)
ATDM_SET_ENABLE(Piro_MatrixFreeDecorator_UnitTests_MPI_4_DISABLE ON)

# Disable Belos resolve tests that are randomly failing (#2965)
ATDM_SET_ENABLE(Belos_resolve_cg_hb_MPI_4_DISABLE ON)
ATDM_SET_ENABLE(Belos_resolve_gmres_hb_1_MPI_4_DISABLE ON)
