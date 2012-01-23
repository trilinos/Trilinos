
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${TEUCHOS_STD_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  NormalsEvaluator
  SOURCES 
     normals.cpp 
     ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  IntegratorScalar 
  SOURCES 
     integrator_scalar.cpp 
     ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  GatherCoordinates 
  SOURCES 
     gather_coordinates.cpp 
     ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DOF_PointFields 
  SOURCES 
     dof_pointfield.cpp 
     ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
