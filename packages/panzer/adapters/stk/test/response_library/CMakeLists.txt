
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${${PROJECT_NAME}_HOME_DIR}/Panzer/panzer/test/equation_set)

INCLUDE_DIRECTORIES(${${PROJECT_NAME}_HOME_DIR}/Panzer/panzer/test/closure_model)
INCLUDE_DIRECTORIES(${${PROJECT_NAME}_HOME_DIR}/Panzer/panzer/test/response_library)

INCLUDE_DIRECTORIES(${${PROJECT_NAME}_HOME_DIR}/Panzer/panzer/adapters/stk/test/bcstrategy)

SET(UNIT_TEST_DRIVER
  ${${PROJECT_NAME}_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  STK_ResponseLibraryTest
  SOURCES response_library.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )
  
