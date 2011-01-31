
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  XML_test  
  SOURCES cxx_main.cpp 
  COMM serial mpi  
  NUM_MPI_PROCS 1
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  XMLParameterListHelpers_UnitTest
  SOURCES XMLParameterListHelpers_UnitTest.cpp
    ../UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
