TRIBITS_ADD_EXECUTABLE(
	frosch_test_import
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
    frosch_test_import
    NAME frosch_test_import
    NUM_MPI_PROCS 8
)