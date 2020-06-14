TRIBITS_ADD_EXECUTABLE(
  Test_InterpolationProjection_HEX
  SOURCES test_interpolation_projection_HEX.cpp
  ADD_DIR_TO_NAME
  )

TRIBITS_ADD_TEST(
  Test_InterpolationProjection_HEX
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  RUN_SERIAL
  )


TRIBITS_ADD_EXECUTABLE(
  Test_InterpolationProjection_TET
  SOURCES test_interpolation_projection_TET.cpp
  ADD_DIR_TO_NAME
  )

TRIBITS_ADD_TEST(
  Test_InterpolationProjection_TET
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  RUN_SERIAL
  )

    
TRIBITS_ADD_EXECUTABLE(
  Test_InterpolationProjection_QUAD
  SOURCES test_interpolation_projection_QUAD.cpp
  ADD_DIR_TO_NAME
  )
  
TRIBITS_ADD_TEST(
  Test_InterpolationProjection_QUAD
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  RUN_SERIAL
  )


TRIBITS_ADD_EXECUTABLE(
  Test_InterpolationProjection_TRI
  SOURCES test_interpolation_projection_TRI.cpp
  ADD_DIR_TO_NAME
  )
  
TRIBITS_ADD_TEST(
  Test_InterpolationProjection_TRI
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  RUN_SERIAL
  )


TRIBITS_ADD_EXECUTABLE(
  Test_Convergence_HEX
  SOURCES test_convergence_HEX.cpp
  ADD_DIR_TO_NAME
  )
  
TRIBITS_ADD_TEST(
  Test_Convergence_HEX
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  RUN_SERIAL 
  )


TRIBITS_ADD_EXECUTABLE(
  Test_Convergence_TET
  SOURCES test_convergence_TET.cpp
  ADD_DIR_TO_NAME
  )
  
TRIBITS_ADD_TEST(
  Test_Convergence_TET
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  RUN_SERIAL 
  )


TRIBITS_ADD_EXECUTABLE(
  Test_Convergence_QUAD
  SOURCES test_convergence_QUAD.cpp
  ADD_DIR_TO_NAME
  )
  
TRIBITS_ADD_TEST(
  Test_Convergence_QUAD
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  RUN_SERIAL
  )


TRIBITS_ADD_EXECUTABLE(
  Test_Convergence_TRI
  SOURCES test_convergence_TRI.cpp
  ADD_DIR_TO_NAME
  )

TRIBITS_ADD_TEST(
  Test_Convergence_TRI
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  RUN_SERIAL
  )
