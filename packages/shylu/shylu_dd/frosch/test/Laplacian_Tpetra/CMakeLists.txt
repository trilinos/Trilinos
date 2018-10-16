TRIBITS_ADD_EXECUTABLE(
    laplacian_tpetra
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(TpetraCopyFiles
    DEST_FILES ParametersGDSW.xml ParametersRGDSW.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS laplacian_tpetra
)

#AH 04/05/2018: Adding Tests
TRIBITS_ADD_TEST(
    laplacian_tpetra
    NAME test_laplacian_tpetra_2d_gdsw
    ARGS "--DIM=2"
    COMM mpi
    NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
    laplacian_tpetra
    NAME test_laplacian_tpetra_3d_gdsw
    ARGS "--DIM=3"
    COMM mpi
    NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
    laplacian_tpetra
    NAME test_laplacian_tpetra_2d_rgdsw
    ARGS "--DIM=2 --RGDSW"
    COMM mpi
    NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
    laplacian_tpetra
    NAME test_laplacian_tpetra_3d_rgdsw
    ARGS "--DIM=3 --RGDSW"
    COMM mpi
    NUM_MPI_PROCS 8
)
