
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SortManager_test
  SOURCES cxx_main.cpp
  COMM serial mpi
  ARGS
    "--verbose"
#    "--debug"
  STANDARD_PASS_OUTPUT
  XHOST KALLIKRATES
  )
