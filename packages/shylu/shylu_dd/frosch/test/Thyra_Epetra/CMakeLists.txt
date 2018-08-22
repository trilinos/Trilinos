TRIBITS_ADD_EXECUTABLE(
    thyra_epetra
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(ThyraEpetraCopyFiles
    DEST_FILES stratimikos_ParameterList.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS thyra_epetra
)

TRIBITS_ADD_TEST(
    thyra_epetra
    NAME test_thyra_epetra
    COMM mpi
    NUM_MPI_PROCS 4
)