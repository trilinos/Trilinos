
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_OrientationEncoding
  SOURCES test_orientation_encoding.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  " 
  ADD_DIR_TO_NAME
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_OrientationToolsCoeffMatrix
  SOURCES test_orientationtools_coeff_matrix.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  " 
  ADD_DIR_TO_NAME
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_OrientationToolsMesh
  SOURCES test_orientationtools_mesh.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  " 
  ADD_DIR_TO_NAME
  )

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

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_Orientation_QUAD
  SOURCES test_orientation_QUAD.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )  

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_Orientation_TRI
  SOURCES test_orientation_TRI.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  ) 
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_Orientation_TET
  SOURCES test_orientation_TET.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )
   
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_Orientation_HEX
  SOURCES test_orientation_HEX.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )
