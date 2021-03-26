include_directories(${ROL2_TYPEU_TESTPROBLEMS})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  trust_region_algorithm
  SOURCES trust_region_algorithm.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
)

