
INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Time_test  
  SOURCES cxx_main.cpp 
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS -v
  STANDARD_PASS_OUTPUT
  )

# The TimeMonitor unit tests, in particular the summarize() test,
# aren't useful unless run with more than one MPI process.  We pick 3
# because it's odd (thereby testing possible corner cases of the
# parallel recursions in the current summarize() implementation) and
# bigger than 1.
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TimeMonitor_UnitTests
  SOURCES
    TimeMonitor_UnitTests.cpp
    ../UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 3
  ARGS --teuchos-suppress-startup-banner
  XHOSTTYPE_TEST AIX
  STANDARD_PASS_OUTPUT
  )
