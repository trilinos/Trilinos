

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/disc-fe/test/equation_set)
TRIBITS_INCLUDE_DIRECTORIES(${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/disc-fe/test/closure_model)
TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/test/bcstrategy)

SET(UNIT_TEST_DRIVER
  ${PANZER_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  initial_condition_builder
  SOURCES initial_condition_builder.cpp user_app_STKClosureModel_Factory.hpp user_app_STKClosureModel_Factory_impl.hpp user_app_STKClosureModel_Factory_TemplateBuilder.hpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  initial_condition_builder2
  SOURCES initial_condition_builder2.cpp user_app_STKClosureModel_Factory.hpp user_app_STKClosureModel_Factory_impl.hpp user_app_STKClosureModel_Factory_TemplateBuilder.hpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  initial_condition_control
  SOURCES initial_condition_control.cpp user_app_STKClosureModel_Factory.hpp user_app_STKClosureModel_Factory_impl.hpp user_app_STKClosureModel_Factory_TemplateBuilder.hpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(ic_builder2_files
  SOURCE_FILES block-decomp.exo block-decomp.exo.2.0 block-decomp.exo.2.1
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/"
  EXEDEPS initial_condition_builder
  )
