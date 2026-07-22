TRIBITS_COPY_FILES_TO_BINARY_DIR(
  Matlab_CustomP_Tests_cp
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCE_FILES
    customP.m
    CustomPTest.m
    laplacianfun.m
    matlabParams.xml
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

TRIBITS_ADD_ADVANCED_TEST(
  Matlab_CustomP
  TEST_0
  CMND ${CMAKE_CURRENT_BINARY_DIR}/../../bin/matlab
  ARGS -r CustomPTest
  PASS_REGULAR_EXPRESSION "passed"
  CATEGORIES BASIC NIGHTLY
  COMM serial
)
