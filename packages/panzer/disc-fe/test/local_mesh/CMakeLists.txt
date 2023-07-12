
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${PANZER_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LocalMesh_Tests
  SOURCES setup_sublocal_mesh_info.cpp 
          local_partitioning_utilities.cpp
          setup_partitioned_workset_utilities.cpp
          ${UNIT_TEST_DRIVER}
  COMM mpi
  NUM_MPI_PROCS 1
  )
