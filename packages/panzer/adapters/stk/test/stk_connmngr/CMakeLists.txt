
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${TEUCHOS_STD_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tSTKConnManager
  SOURCES tSTKConnManager.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tSquareQuadMeshDOFManager
  SOURCES tSquareQuadMeshDOFManager.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tSquareQuadMeshDOFManager_FEI
  SOURCES tSquareQuadMeshDOFManager_FEI.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tDOFManager2_Orientation
  SOURCES tDOFManager2_Orientation.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 2
  COMM serial mpi
  )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tSquareTriMeshDOFManager
  SOURCES tSquareTriMeshDOFManager.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tSquareTriMeshDOFManager_FEI
  SOURCES tSquareTriMeshDOFManager_FEI.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tEpetraLinObjFactory
  SOURCES tEpetraLinearObjFactory.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tCubeHexMeshDOFManager
  SOURCES tCubeHexMeshDOFManager.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tSquareQuadMeshDOFManager_edgetests
  SOURCES tSquareQuadMeshDOFManager_edgetests.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tSquareQuadMeshDOFManager_edgetests_FEI
  SOURCES tSquareQuadMeshDOFManager_edgetests_FEI.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tBlockedDOFManagerFactory
  SOURCES tBlockedDOFManagerFactory.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tDOFManager2_SimpleTests
  SOURCES tDOFManager2_SimpleTests.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 4
  COMM mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SingleBlockDOFManagerTests
  SOURCES SingleBlockDOFManagerTests.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 3
  COMM mpi
  )
