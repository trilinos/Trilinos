
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp)

SET(SOURCES
    equation_set.cpp 
    user_app_EquationSetFactory.hpp 
    user_app_EquationSet_Energy.hpp
    user_app_EquationSet_EnergyT.hpp
    ) 

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  equation_set
  SOURCES ${SOURCES} ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
