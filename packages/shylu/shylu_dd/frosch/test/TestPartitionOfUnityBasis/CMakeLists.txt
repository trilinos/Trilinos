TRIBITS_ADD_EXECUTABLE(
	os_partitionofunitybasis
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	os_partitionofunitybasis
	NAME test_os_partitionofunitybasis
    COMM mpi
    NUM_MPI_PROCS 1
)
