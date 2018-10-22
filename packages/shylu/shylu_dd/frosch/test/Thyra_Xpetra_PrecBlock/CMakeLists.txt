TRIBITS_ADD_EXECUTABLE(
    thyra_xpetra_precblock
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(ThyraXpetraPrecBlockCopyFiles
    DEST_FILES xpetra_ParameterList.xml stokes.h5
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS thyra_xpetra_prec
)

TRIBITS_ADD_TEST(
    thyra_xpetra_precblock
    NAME test_thyra_xpetra_precblock
    COMM mpi
    NUM_MPI_PROCS 4
)
