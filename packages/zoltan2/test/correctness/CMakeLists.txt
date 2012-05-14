#
# RCB tests done by Zoltan in nightly testing
#

INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  zoltanCompare
  SOURCES zoltanCompare.cpp
  NUM_MPI_PROCS 2
  COMM mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )
