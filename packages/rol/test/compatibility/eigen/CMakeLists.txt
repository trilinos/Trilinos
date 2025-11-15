# Eigen compatibility tests

ROL_ADD_EXECUTABLE_AND_TEST(
  eigen_LapackCompatibility
  SOURCES test_01.cpp
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
)