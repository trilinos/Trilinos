TRIBITS_ADD_EXECUTABLE(
    thyra_xpetra
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(ThyraXpetraCopyFiles
    DEST_FILES xpetra_ParameterList.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS thyra_xpetra
)

TRIBITS_ADD_TEST(
    thyra_xpetra
    NAME test_thyra_xpetra
    COMM mpi
    NUM_MPI_PROCS 4
)