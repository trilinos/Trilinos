
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_OrientationToolsModifyBasis_LowerOrder
  SOURCES test_orientationtools_modify_basis_lower_order.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  " 
  ADD_DIR_TO_NAME
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_OrientationToolsModifyBasis_HGRAD
  SOURCES test_orientationtools_modify_basis_hgrad.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  " 
  ADD_DIR_TO_NAME
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_OrientationToolsModifyBasis_HCURL
  SOURCES test_orientationtools_modify_basis_hcurl.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  " 
  ADD_DIR_TO_NAME
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_OrientationToolsModifyBasis_HDIV
  SOURCES test_orientationtools_modify_basis_hdiv.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  " 
  ADD_DIR_TO_NAME
  )

