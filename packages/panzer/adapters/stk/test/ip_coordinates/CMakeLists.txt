
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${${PROJECT_NAME}_SOURCE_DIR}/Panzer/panzer/test/equation_set)

INCLUDE_DIRECTORIES(${${PROJECT_NAME}_SOURCE_DIR}/Panzer/panzer/test/closure_model)
INCLUDE_DIRECTORIES(${${PROJECT_NAME}_SOURCE_DIR}/Panzer/panzer/test/response_library)

INCLUDE_DIRECTORIES(${${PROJECT_NAME}_SOURCE_DIR}/Panzer/panzer/adapters/stk/test/bcstrategy)

SET(UNIT_TEST_DRIVER
  ${TEUCHOS_STD_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ip_coordinates
  SOURCES ip_coordinates.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )
  
