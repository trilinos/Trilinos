TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${PANZER_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tCartesianDOFMgr_DynRankView
  SOURCES tCartesianDOFMgr_DynRankView.cpp CartesianConnManager.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 4
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tCartesianDOFMgr_HighOrder
  SOURCES tCartesianDOFMgr_HighOrder.cpp CartesianConnManager.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 4
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tCartesianDOFMgr_DG
  SOURCES tCartesianDOFMgr_DG.cpp CartesianConnManager.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 4
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tCartesianTop
  SOURCES tCartesianTop.cpp CartesianConnManager.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 4
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tCartesianDOFMgr
  SOURCES tCartesianDOFMgr.cpp CartesianConnManager.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 4
  COMM serial mpi
  )
