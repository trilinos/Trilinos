TRIBITS_ADD_EXECUTABLE(
	interfacepartitionofunity
	SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(IPOUCopyFiles
    DEST_FILES ParametersIPOU.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS interfacepartitionofunity
)

IF(ShyLU_DDFROSch_ENABLE_Epetra)
TRIBITS_ADD_TEST(
interfacepartitionofunity
NAME test_interfacepartitionofunity_DIM2_EPETRA
ARGS "--DIM=2 --M=4 --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
interfacepartitionofunity
NAME test_interfacepartitionofunity_DIM3_EPETRA
ARGS "--DIM=3 --M=4 --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

TRIBITS_ADD_TEST(
interfacepartitionofunity
NAME test_interfacepartitionofunity_DIM2_TPETRA
ARGS "--DIM=2 --M=4 --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
interfacepartitionofunity
NAME test_interfacepartitionofunity_DIM3_TPETRA
ARGS "--DIM=3 --M=4 --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
