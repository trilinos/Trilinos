
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  periodic_bcs 
  SOURCES periodic_bcs.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 4
  )
  
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  periodic_mesh
  SOURCES periodic_mesh.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )
