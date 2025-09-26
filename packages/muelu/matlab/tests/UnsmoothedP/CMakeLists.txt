TRIBITS_COPY_FILES_TO_BINARY_DIR(
  Matlab_Ptent_Tests_cp
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCE_FILES
    createP.m
    laplacianfun.m
    matlabParams.xml
    matlabParamsEasy.xml
    mueluParams.xml
    Ptent_UnitTest.m
    PtentEasy_UnitTest.m
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

TRIBITS_ADD_ADVANCED_TEST(
  Matlab_Ptent
  TEST_0
  CMND ${CMAKE_CURRENT_BINARY_DIR}/../../bin/matlab
  ARGS -r Ptent_UnitTest
  PASS_REGULAR_EXPRESSION "passed"
  CATEGORIES BASIC NIGHTLY
  COMM serial
)

TRIBITS_ADD_ADVANCED_TEST(
  Matlab_PtentEasy
  TEST_0
  CMND ${CMAKE_CURRENT_BINARY_DIR}/../../bin/matlab
  ARGS -r PtentEasy_UnitTest
  PASS_REGULAR_EXPRESSION "passed"
  CATEGORIES BASIC NIGHTLY
  COMM serial
)
