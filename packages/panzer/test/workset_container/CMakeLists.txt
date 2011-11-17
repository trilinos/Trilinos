
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${${PROJECT_NAME}_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp)

SET(SOURCES
    workset_container.cpp 
    ) 

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  workset_container 
  SOURCES ${SOURCES} ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
