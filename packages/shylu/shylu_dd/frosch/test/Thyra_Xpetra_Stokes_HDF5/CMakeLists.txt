TRIBITS_ADD_EXECUTABLE(
    thyra_xpetra_stokes_hdf5
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(ThyraXpetraPrecBlockCopyFiles
    DEST_FILES ParameterList_TwoLevelBlockPreconditioner_GDSW.xml stokes.h5
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS thyra_xpetra_stokes_hdf5
)

####################################
# TwoLevelBlockPreconditioner_GDSW #
####################################
## Tpetra
TRIBITS_ADD_TEST(
thyra_xpetra_stokes_hdf5
NAME test_thyra_xpetra_stokes_hdf5_TLBP_GDSW_O0_TPETRA
ARGS "--O=0 --PLIST=ParameterList_TwoLevelBlockPreconditioner_GDSW.xml"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
thyra_xpetra_stokes_hdf5
NAME test_thyra_xpetra_stokes_hdf5_TLBP_GDSW_O1_TPETRA
ARGS "--O=1 --PLIST=ParameterList_TwoLevelBlockPreconditioner_GDSW.xml"
COMM mpi
NUM_MPI_PROCS 4
)
