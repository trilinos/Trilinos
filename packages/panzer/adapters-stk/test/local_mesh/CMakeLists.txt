

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${PANZER_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LocalMesh_Tests
  SOURCES local_mesh_utilities.cpp 
          PanzerSTK_UnitTest_BuildMesh.cpp
          PanzerSTK_UnitTest_STKInterfaceGenerator.cpp
          ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
