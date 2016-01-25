## This test compares Zoltan and Zoltan2
## Don't run this test from tarball release, as not all Zoltan input files
## are available in tarball.

IF (${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    zoltanCompare
    SOURCES zoltanCompare.cpp
    NUM_MPI_PROCS 8 
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()
