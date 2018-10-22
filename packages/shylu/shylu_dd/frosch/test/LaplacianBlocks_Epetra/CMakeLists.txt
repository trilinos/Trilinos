TRIBITS_ADD_EXECUTABLE(
    frosch_laplacianblocks_epetra
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(ShyLU_DDFROSchEpetraBlocksCopyFiles
    DEST_FILES ParametersGDSW.xml ParametersRGDSW.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS frosch_laplacianblocks_epetra
)

TRIBITS_ADD_TEST(
    frosch_laplacianblocks_epetra
    NAME test_frosch_laplacianblocks_epetra_2d_gdsw
    ARGS "--DIM=2"
    COMM mpi
    NUM_MPI_PROCS 4
)
