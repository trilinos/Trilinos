TRIBITS_COPY_FILES_TO_BINARY_DIR(Matlab_Reentrant_Tests_cp
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES modifyP.m laplacianfun.m matlabParams.xml Reentrant_Test.m
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    )

TRIBITS_ADD_ADVANCED_TEST(
  Matlab_Reentrant
  TEST_0 CMND ${MATLAB_ROOT}/bin/matlab ARGS -nodesktop -nosplash -r Reentrant_Test
  PASS_REGULAR_EXPRESSION "passed"
  CATEGORIES BASIC NIGHTLY
  COMM serial
  )
