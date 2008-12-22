
INCLUDE(Trilinos_Add_Executable_And_Test)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  performance_bracket_operator
  SOURCES Performance_BracketOperator.cpp
  DEPLIBS phalanx_test_utilities
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Test passed!"
  )