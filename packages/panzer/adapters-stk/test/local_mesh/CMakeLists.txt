

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

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

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_mesh.exo.2.0 ${CMAKE_CURRENT_BINARY_DIR}/test_mesh.exo.2.0 COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_mesh.exo.2.1 ${CMAKE_CURRENT_BINARY_DIR}/test_mesh.exo.2.1 COPYONLY)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LocalMesh_ParallelTests
  SOURCES parallel_local_mesh_utilities.cpp
          PanzerSTK_UnitTest_BuildMesh.cpp
          PanzerSTK_UnitTest_STKInterfaceGenerator.cpp 
          ${UNIT_TEST_DRIVER}
  COMM mpi
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LocalMesh_ParallelPeriodicTests
  SOURCES parallel_periodic_local_mesh_utilities.cpp 
          PanzerSTK_UnitTest_BuildMesh.cpp
          PanzerSTK_UnitTest_STKInterfaceGenerator.cpp
          ${UNIT_TEST_DRIVER}
  COMM mpi
  NUM_MPI_PROCS 2
  )
  