TRIBITS_ADD_EXECUTABLE(
	os_repeatedmap
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	os_repeatedmap
	NAME test_os_repeatedmap_2D
	ARGS "--DIM=2 --M=4"
	COMM mpi
	NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
    os_repeatedmap
    NAME test_os_repeatedmap_3D
    ARGS "--DIM=3 --M=4"
    COMM mpi
    NUM_MPI_PROCS 9
)