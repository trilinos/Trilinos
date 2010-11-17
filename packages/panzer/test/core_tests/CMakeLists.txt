
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  integration_rule
  SOURCES integration_rule.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  integration_values
  SOURCES integration_values.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
  
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  dimension
  SOURCES dimension.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  version
  SOURCES version.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
  
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  basis
  SOURCES basis.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
  
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  basis_values
  SOURCES basis_values.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
      
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  input_equation_set
  SOURCES input_equation_set.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
   
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  input_physics_block
  SOURCES input_physics_block.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
 
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  boundary_condition
  SOURCES boundary_condition.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
  