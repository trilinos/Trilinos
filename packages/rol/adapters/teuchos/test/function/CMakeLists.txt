
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TeuchosFunctionInterface
  SOURCES test_01.cpp
  ARGS PrintItAll
  COMM serial mpi  
  PASS_REGULAR_EXPRESSION "TEST PASSED"
)

