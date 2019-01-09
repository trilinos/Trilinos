TRIBITS_ADD_EXECUTABLE(
	interfacesets
	SOURCES main.cpp
)

IF(ShyLU_DDFROSch_ENABLE_Epetra)
TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM2_EPETRA
ARGS "--DIM=2 --M=4 --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM3_EPETRA
ARGS "--DIM=3 --M=4 --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM2_TPETRA
ARGS "--DIM=2 --M=4 --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM3_TPETRA
ARGS "--DIM=3 --M=4 --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
