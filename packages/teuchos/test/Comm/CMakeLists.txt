
# Required modules
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  DefaultMpiComm_UnitTests  
  SOURCES
    DefaultMpiComm_UnitTests.cpp
    ../UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE(
  ReduceAllAndScatter_test  
  COMM serial mpi  
  SOURCES ReduceAllAndScatter_test.cpp
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Comm_test  
  SOURCES Comm_test.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1-10
  STANDARD_PASS_OUTPUT
  )
