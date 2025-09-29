TRIBITS_COPY_FILES_TO_BINARY_DIR(
  Matlab_Jacobi_Test_cp
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCE_FILES
    getDiag.m
    jacobi.m
    JacobiSmoother_UnitTest.m
    laplacianfun.m
    matlabParams.xml
    mueluParams.xml
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

TRIBITS_ADD_ADVANCED_TEST(
  Matlab_Jacobi
  TEST_0
  CMND ${CMAKE_CURRENT_BINARY_DIR}/../../bin/matlab
  ARGS -r JacobiSmoother_UnitTest
  PASS_REGULAR_EXPRESSION "passed"
  CATEGORIES BASIC NIGHTLY
  COMM serial
)
