
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsAddAdvancedTest)


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CWrapperSupport_UnitTests
  SOURCES
    Teuchos_CWrapperSupport_UnitTests.cpp
    someCFunc.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
