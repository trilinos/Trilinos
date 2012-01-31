
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${${PROJECT_NAME}_SOURCE_DIR}/Panzer/panzer/test/equation_set)

INCLUDE_DIRECTORIES(${${PROJECT_NAME}_SOURCE_DIR}/Panzer/panzer/test/closure_model)

INCLUDE_DIRECTORIES(${${PROJECT_NAME}_SOURCE_DIR}/Panzer/panzer/adapters/stk/test/bcstrategy)

SET(UNIT_TEST_DRIVER
  ${TEUCHOS_STD_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  initial_condition_builder
  SOURCES initial_condition_builder.cpp user_app_STKClosureModel_Factory.hpp user_app_STKClosureModel_Factory_impl.hpp user_app_STKClosureModel_Factory_TemplateBuilder.hpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
  
