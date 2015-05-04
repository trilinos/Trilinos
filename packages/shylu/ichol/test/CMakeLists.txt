
TRIBITS_ADD_EXECUTABLE_AND_TEST(
    test_serial
    SOURCES main_serial.cpp 
    COMM serial mpi
    NUM_MPI_PROCS 1
    FAIL_REGULAR_EXPRESSION "FAILED"
)

