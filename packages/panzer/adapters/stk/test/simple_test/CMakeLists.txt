
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${${PROJECT_NAME}_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tVersion
  SOURCES tVersion.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  )
