include(TribitsAddExecutable)
include(TribitsAddAdvancedTest)

TRIBITS_ADD_EXECUTABLE(
  a_test
  SOURCES a_test.cpp
)

TRIBITS_ADD_ADVANCED_TEST(
  test_of_a
  TEST_0 EXEC a_test
  PASS_REGULAR_EXPRESSION_ALL 
      "A label is: A"
      "A deps are:"
)
