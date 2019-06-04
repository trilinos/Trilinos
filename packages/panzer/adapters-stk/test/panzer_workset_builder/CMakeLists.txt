

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PARENT_PACKAGE_SOURCE_DIR}/disc-fe/test/equation_set)

SET(UNIT_TEST_DRIVER
  ${PANZER_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  workset_builder
  SOURCES workset_builder.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  d_workset_builder
  SOURCES d_workset_builder.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  d_workset_builder_3d
  SOURCES d_workset_builder_3d.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  cascade 
  SOURCES cascade.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  hdiv_basis 
  SOURCES hdiv_basis.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  workset_container 
  SOURCES workset_container.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  check_rotation_matrices 
  SOURCES check_rotation_matrices.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  check_sorted_face_quads
  SOURCES check_sorted_face_quads.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
