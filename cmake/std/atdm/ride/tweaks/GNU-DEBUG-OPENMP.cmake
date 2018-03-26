# #2410: This passes in the GNU-RELESE-OPENMP build so it is okay to disable
# in this DEBUG build.  This is likey related to the mix-lanauage compiler
# defect reported in #1208.
ATDM_SET_ENABLE(TeuchosNumerics_LAPACK_test_MPI_1_DISABLE ON)

# Disable test that runs over 30 min currently (#2446)
ATDM_SET_ENABLE(PanzerAdaptersSTK_MixedPoissonExample-ConvTest-Hex-Order-3_DISABLE ON)
