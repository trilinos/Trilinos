TRIBITS_ADD_EXECUTABLE(
	interfacesets
	SOURCES main.cpp
)

############################
# Communication Strategy 0 #
############################
TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM2_CSCommCrsMatrix_TPETRA
ARGS "--DIM=2 --M=4 --COMMSTRAT=0"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM3_CSCommCrsMatrix_TPETRA
ARGS "--DIM=3 --M=4 --COMMSTRAT=0"
COMM mpi
NUM_MPI_PROCS 8
)

############################
# Communication Strategy 1 #
############################
TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM2_CSCommCrsGraph_TPETRA
ARGS "--DIM=2 --M=4 --COMMSTRAT=1"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM3_CSCommCrsGraph_TPETRA
ARGS "--DIM=3 --M=4 --COMMSTRAT=1"
COMM mpi
NUM_MPI_PROCS 8
)

############################
# Communication Strategy 2 #
############################
TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM2_CSOneToOne_TPETRA
ARGS "--DIM=2 --M=4 --COMMSTRAT=2"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM3_CSOneToOne_TPETRA
ARGS "--DIM=3 --M=4 --COMMSTRAT=2"
COMM mpi
NUM_MPI_PROCS 8
)
