INCLUDE(TribitsAddExecutableAndTest)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FunctionObject_Tests
  SOURCES cxx_main.cpp FunctionObject_UnitTests.cpp
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FunctionObject_Serialization_Tests
  SOURCES cxx_main.cpp FunctionObject_Serialization_UnitTests.cpp
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )



