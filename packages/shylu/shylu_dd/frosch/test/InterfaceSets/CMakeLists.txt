TRIBITS_ADD_EXECUTABLE(
	interfacesets
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	interfacesets
	NAME test_interfacesets_2D
	ARGS "--DIM=2 --M=4"
	COMM mpi
	NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
    repeatedmap
    NAME test_interfacesets_3D
    ARGS "--DIM=3 --M=4"
    COMM mpi
    NUM_MPI_PROCS 8
)
