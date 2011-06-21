
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/equation_set)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/closure_model)

SET(UNIT_TEST_DRIVER
  ${TRILINOS_HOME_DIR}/packages/teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp)

SET(SOURCES
    bcstrategy.cpp
    user_app_BCStrategy_Dirichlet_Constant.hpp
    user_app_BCStrategy_Dirichlet_ConstantT.hpp
    user_app_BCStrategy_Factory.hpp
    ) 

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  bcstrategy
  SOURCES ${SOURCES} ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
