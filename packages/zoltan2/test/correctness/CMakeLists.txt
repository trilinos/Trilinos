#
# Tests that check correctness of RCB.
#

INCLUDE(TribitsAddExecutableAndTest)


## This test calls zdrive functions to read
## in the same input data and run some of
## the same tests that Zoltan runs.

IF (HAVE_${PACKAGE_NAME}_ZOLTAN)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  zoltanCompare
  SOURCES zoltanCompare.cpp
  NUM_MPI_PROCS 2 3 4 5 6 8
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

ENDIF()
