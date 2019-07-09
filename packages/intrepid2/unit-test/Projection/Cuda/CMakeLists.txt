TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_InterpolationProjection_HEX
  SOURCES test_interpolation_projection_HEX.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_InterpolationProjection_TET
  SOURCES test_interpolation_projection_TET.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )
    
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_InterpolationProjection_QUAD
  SOURCES test_interpolation_projection_QUAD.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_InterpolationProjection_TRI
  SOURCES test_interpolation_projection_TRI.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )