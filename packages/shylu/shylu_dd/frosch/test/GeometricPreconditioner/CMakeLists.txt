TRIBITS_ADD_EXECUTABLE(
	geometricpreconditioner
	SOURCES main.cpp
)

###########################
# GeometricPreconditioner #
###########################
## Tpetra
TRIBITS_ADD_TEST(
geometricpreconditioner
NAME test_geometricpreconditioner
ARGS ""
COMM mpi
NUM_MPI_PROCS 2
)
