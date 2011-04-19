
INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Time_test  
  SOURCES cxx_main.cpp 
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS -v
  STANDARD_PASS_OUTPUT
  )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TimeMonitor_UnitTests
  SOURCES
    TimeMonitor_UnitTests.cpp
    ../UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS --teuchos-suppress-startup-banner
  XHOSTTYPE_TEST AIX
  STANDARD_PASS_OUTPUT
  )
