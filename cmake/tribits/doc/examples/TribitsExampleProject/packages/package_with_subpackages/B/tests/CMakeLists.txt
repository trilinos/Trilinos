include(TribitsAddExecutable)
include(TribitsAddAdvancedTest)

TRIBITS_ADD_EXECUTABLE(
  b_test
  SOURCES b_test.cpp
)

TRIBITS_ADD_ADVANCED_TEST(
  test_of_b
  TEST_0 EXEC b_test
  PASS_REGULAR_EXPRESSION_ALL 
      "A label is: A"
      "B label is: B"
      "B deps are: A"
)
