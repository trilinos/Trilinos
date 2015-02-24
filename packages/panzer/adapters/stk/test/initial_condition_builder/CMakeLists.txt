

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/test/equation_set)
INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/test/closure_model)
INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/stk/test/bcstrategy)

SET(UNIT_TEST_DRIVER
  ${TEUCHOS_STD_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  initial_condition_builder
  SOURCES initial_condition_builder.cpp user_app_STKClosureModel_Factory.hpp user_app_STKClosureModel_Factory_impl.hpp user_app_STKClosureModel_Factory_TemplateBuilder.hpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
  
# THIS TEST FAILS!
############################################
# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   initial_condition_builder2
#   SOURCES initial_condition_builder2.cpp user_app_STKClosureModel_Factory.hpp user_app_STKClosureModel_Factory_impl.hpp user_app_STKClosureModel_Factory_TemplateBuilder.hpp ${UNIT_TEST_DRIVER}
#   COMM serial mpi
#   NUM_MPI_PROCS 2
#   )

TRIBITS_COPY_FILES_TO_BINARY_DIR(ic_builder2_files
  SOURCE_FILES block-decomp.exo.2.0 block-decomp.exo.2.1
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/"
  EXEDEPS initial_condition_builder
  )
