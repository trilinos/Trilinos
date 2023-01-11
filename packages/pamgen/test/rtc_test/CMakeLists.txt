
TRIBITS_INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/rtcompiler)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  rtc_test
  SOURCES rtcompiler-test.cc
  COMM serial
  NUM_MPI_PROCS 1
  )

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/table.dat table.dat COPYONLY)
