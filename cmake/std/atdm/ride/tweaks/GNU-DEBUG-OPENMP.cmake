# #2410: This passes in the GNU-RELESE-OPENMP build so it is okay to disable
# in this DEBUG build.  This is likey related to the mix-lanauage compiler
# defect reported in #1208.
ATDM_SET_ENABLE(TeuchosNumerics_LAPACK_test_MPI_1_DISABLE ON)
