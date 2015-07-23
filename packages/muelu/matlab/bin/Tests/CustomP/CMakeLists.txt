 TRIBITS_COPY_FILES_TO_BINARY_DIR(Matlab_CustomP_Tests_cp
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES customP.m laplacianfun.m matlabParams.xml CustomPTest.m
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    )

TRIBITS_ADD_ADVANCED_TEST(
  Matlab_CustomP
  TEST_0 CMND ${MATLAB_ROOT}/bin/matlab ARGS -nodesktop -nosplash -r CustomPTest
  PASS_REGULAR_EXPRESSION "passed"
  CATEGORIES BASIC NIGHTLY
  COMM serial
  )
