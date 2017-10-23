TRIBITS_ADD_EXECUTABLE(
    frosch_laplacian_epetra
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(ShyLU_DDFROSchEpetraCopyFiles
    DEST_FILES ParametersGDSW.xml ParametersRGDSW.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS frosch_laplacian_epetra
)

TRIBITS_ADD_TEST(
    frosch_laplacian_epetra
    NAME test_frosch_laplacian_epetra
    COMM mpi
    NUM_MPI_PROCS 9
)