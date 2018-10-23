TRIBITS_ADD_EXECUTABLE(
    thyra_tpetra
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(ThyraTpetraCopyFiles
    DEST_FILES xpetra_ParameterList.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS thyra_tpetra
)

TRIBITS_ADD_TEST(
    thyra_tpetra
    NAME test_thyra_tpetra
    COMM mpi
    NUM_MPI_PROCS 4
)