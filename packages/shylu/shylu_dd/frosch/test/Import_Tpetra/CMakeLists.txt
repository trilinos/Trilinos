TRIBITS_ADD_EXECUTABLE(
	test_import
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
    test_import
    NAME test_import
    NUM_MPI_PROCS 8
)
