
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  XML_test  
  SOURCES cxx_main.cpp 
  COMM serial mpi  
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  XMLParameterListHelpers_UnitTest
  SOURCES 
    XMLParameterListHelpers_UnitTest.cpp
    XMLParameterListHelpers_AnonSublists.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
