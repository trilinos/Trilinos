TRIBITS_ADD_EXECUTABLE(
    laplacian_epetra
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(EpetraCopyFiles
    DEST_FILES ParametersGDSW.xml ParametersRGDSW.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS laplacian_epetra
)

#AH 04/05/2018: Adding Tests
TRIBITS_ADD_TEST(
    laplacian_epetra
    NAME test_laplacian_epetra_2d_gdsw
    ARGS "--DIM=2"
    COMM mpi
    NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
    laplacian_epetra
    NAME test_laplacian_epetra_3d_gdsw
    ARGS "--DIM=3"
    COMM mpi
    NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
    laplacian_epetra
    NAME test_laplacian_epetra_2d_rgdsw
    ARGS "--DIM=2 --RGDSW"
    COMM mpi
    NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
    laplacian_epetra
    NAME test_laplacian_epetra_3d_rgdsw
    ARGS "--DIM=3 --RGDSW"
    COMM mpi
    NUM_MPI_PROCS 8
)
