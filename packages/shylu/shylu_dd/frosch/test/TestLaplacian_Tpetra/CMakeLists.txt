TRIBITS_ADD_EXECUTABLE(
    frosch_laplacian_tpetra
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(ShyLU_DDFROSchTpetraCopyFiles
    DEST_FILES ParametersGDSW.xml ParametersRGDSW.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS frosch_laplacian_tpetra
)

#AH 09/20/2017: Removed due to issue in Tpetra
#TRIBITS_ADD_TEST(
#    frosch_laplacian_tpetra
#    NAME test_frosch_laplacian_tpetra
#    COMM mpi
#    NUM_MPI_PROCS 9
#)