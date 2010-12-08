
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${TRILINOS_HOME_DIR}/Panzer/panzer/test/equation_set)

INCLUDE_DIRECTORIES(${TRILINOS_HOME_DIR}/Panzer/panzer/test/model_factory)

INCLUDE_DIRECTORIES(${TRILINOS_HOME_DIR}/Panzer/panzer/test/bcstrategy)

SET(UNIT_TEST_DRIVER
  ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  assembly_engine
  SOURCES assembly_engine.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 4
  )
  