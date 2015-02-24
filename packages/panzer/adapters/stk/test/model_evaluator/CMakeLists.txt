

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/equation_set)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/closure_model)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/adapters/stk/test/bcstrategy)

SET(UNIT_TEST_DRIVER
  ${TEUCHOS_STD_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  model_evaluator
  SOURCES model_evaluator.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 4
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  model_evaluator_mass_check
  SOURCES model_evaluator_mass_check.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  thyra_model_evaluator
  SOURCES thyra_model_evaluator.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 4
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  explicit_model_evaluator
  SOURCES explicit_model_evaluator.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 4
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  response_residual 
  SOURCES response_residual.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )
