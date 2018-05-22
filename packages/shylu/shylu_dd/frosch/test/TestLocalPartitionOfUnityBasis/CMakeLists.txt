TRIBITS_ADD_EXECUTABLE(
	frosch_localpartitionofunitybasis
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	frosch_localpartitionofunitybasis
	NAME test_frosch_localpartitionofunitybasis
    COMM mpi
    NUM_MPI_PROCS 1
)
