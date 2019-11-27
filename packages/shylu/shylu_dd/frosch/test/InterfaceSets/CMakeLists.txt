TRIBITS_ADD_EXECUTABLE(
	interfacesets
	SOURCES main.cpp
)

IF(HAVE_SHYLU_DDFROSCH_EPETRA AND NOT Tpetra_DefaultNode STREQUAL "Kokkos::Compat::KokkosCudaWrapperNode")
############################
# Communication Strategy 0 #
############################
TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM2_CSCommCrsMatrix_EPETRA
ARGS "--DIM=2 --M=4 --COMMSTRAT=0 --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM3_CSCommCrsMatrix_EPETRA
ARGS "--DIM=3 --M=4 --COMMSTRAT=0 --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

############################
# Communication Strategy 1 #
############################
TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM2_CSCommCrsGraph_EPETRA
ARGS "--DIM=2 --M=4 --COMMSTRAT=1 --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM3_CSCommCrsGraph_EPETRA
ARGS "--DIM=3 --M=4 --COMMSTRAT=1 --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

############################
# Communication Strategy 2 #
############################
TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM2_CSOneToOne_EPETRA
ARGS "--DIM=2 --M=4 --COMMSTRAT=2 --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM3_CSOneToOne_EPETRA
ARGS "--DIM=3 --M=4 --COMMSTRAT=2 --USEEPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
ENDIF()

############################
# Communication Strategy 0 #
############################
TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM2_CSCommCrsMatrix_TPETRA
ARGS "--DIM=2 --M=4 --COMMSTRAT=0 --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM3_CSCommCrsMatrix_TPETRA
ARGS "--DIM=3 --M=4 --COMMSTRAT=0 --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

############################
# Communication Strategy 1 #
############################
TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM2_CSCommCrsGraph_TPETRA
ARGS "--DIM=2 --M=4 --COMMSTRAT=1 --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM3_CSCommCrsGraph_TPETRA
ARGS "--DIM=3 --M=4 --COMMSTRAT=1 --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)

############################
# Communication Strategy 2 #
############################
TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM2_CSOneToOne_TPETRA
ARGS "--DIM=2 --M=4 --COMMSTRAT=2 --USETPETRA"
COMM mpi
NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
interfacesets
NAME test_interfacesets_DIM3_CSOneToOne_TPETRA
ARGS "--DIM=3 --M=4 --COMMSTRAT=2 --USETPETRA"
COMM mpi
NUM_MPI_PROCS 8
)
