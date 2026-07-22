TRIBITS_COPY_FILES_TO_BINARY_DIR(
  Matlab_BrickAgg_Tests_cp
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCE_FILES
    BrickAggregationUnitTest.m
    BrickAggregationUnitTestEasy.m
    laplacianfun.m
    matlabParams.xml
    matlabParamsEasy.xml
    mueluParams.xml
    simpleAggregation.m
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

TRIBITS_ADD_ADVANCED_TEST(
  Matlab_BrickAgg
  TEST_0
  CMND ${CMAKE_CURRENT_BINARY_DIR}/../../bin/matlab
  ARGS -r BrickAggregationUnitTest
  PASS_REGULAR_EXPRESSION "passed"
  CATEGORIES BASIC NIGHTLY
  COMM serial
)

TRIBITS_ADD_ADVANCED_TEST(
  Matlab_BrickAggEasy
  TEST_0
  CMND ${CMAKE_CURRENT_BINARY_DIR}/../../bin/matlab
  ARGS -r BrickAggregationUnitTestEasy
  PASS_REGULAR_EXPRESSION "passed"
  CATEGORIES BASIC NIGHTLY
  COMM serial
)
