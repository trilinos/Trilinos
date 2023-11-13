

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_INCLUDE_DIRECTORIES(${PARENT_PACKAGE_SOURCE_DIR}/disc-fe/test/equation_set)
TRIBITS_INCLUDE_DIRECTORIES(${PARENT_PACKAGE_SOURCE_DIR}/disc-fe/test/closure_model)
TRIBITS_INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/bcstrategy)

SET(UNIT_TEST_DRIVER ${PANZER_UNIT_TEST_MAIN})

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

IF(${PARENT_PACKAGE_NAME}_BUILD_HESSIAN_SUPPORT)
   TRIBITS_ADD_EXECUTABLE_AND_TEST(
     model_evaluator_hessians
     SOURCES model_evaluator_hessians.cpp ${UNIT_TEST_DRIVER}
     COMM serial mpi
     NUM_MPI_PROCS 4
     )

   TRIBITS_ADD_EXECUTABLE_AND_TEST(
     model_evaluator_blocked_hessians
     SOURCES model_evaluator_blocked_hessians.cpp ${UNIT_TEST_DRIVER}
     COMM serial mpi
     NUM_MPI_PROCS 4
     )
ENDIF()
