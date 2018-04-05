TRIBITS_ADD_EXECUTABLE(
	frosch_interfacesets
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	frosch_interfacesets
	NAME test_frosch_interfacesets_2D
	ARGS "--DIM=2 --M=4"
	COMM mpi
	NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
    frosch_repeatedmap
    NAME test_frosch_interfacesets_3D
    ARGS "--DIM=3 --M=4"
    COMM mpi
    NUM_MPI_PROCS 8
)
