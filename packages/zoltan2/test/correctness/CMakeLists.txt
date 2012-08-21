#
# Tests that check correctness of RCB.
#

INCLUDE(TribitsAddExecutableAndTest)


## This test calls zdrive functions to read
## in the same input data and run some of
## the same tests that Zoltan runs.

IF (${PACKAGE_NAME}_ENABLE_Experimental)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    zoltanCompare
    SOURCES zoltanCompare.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()
