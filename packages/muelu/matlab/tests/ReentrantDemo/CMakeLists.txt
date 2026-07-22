TRIBITS_COPY_FILES_TO_BINARY_DIR(
  Matlab_Reentrant_Tests_cp
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCE_FILES
    laplacianfun.m
    matlabParams.xml
    modifyP.m
    Reentrant_Test.m
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

TRIBITS_ADD_ADVANCED_TEST(
  Matlab_Reentrant
  TEST_0
  CMND ${CMAKE_CURRENT_BINARY_DIR}/../../bin/matlab
  ARGS -r Reentrant_Test
  PASS_REGULAR_EXPRESSION "passed"
  CATEGORIES BASIC NIGHTLY
  COMM serial
)
