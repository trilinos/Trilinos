
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  tSTKConnManager
  SOURCES tSTKConnManager.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  tSquareQuadMeshDOFManager
  SOURCES tSquareQuadMeshDOFManager.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  tEpetraLinObjFactory
  SOURCES tEpetraLinearObjFactory.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 2
  COMM serial mpi
  )
