
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/equation_set)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/closure_model)

SET(UNIT_TEST_DRIVER
  ${TEUCHOS_STD_UNIT_TEST_MAIN})

SET(SOURCES
    bcstrategy.cpp
    user_app_BCStrategy_Dirichlet_Constant.hpp
    user_app_BCStrategy_Dirichlet_Constant_impl.hpp
    user_app_BCStrategy_Neumann_Constant.hpp
    user_app_BCStrategy_Neumann_Constant_impl.hpp
    user_app_BCStrategy_Factory.hpp
    ) 

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  bcstrategy
  SOURCES ${SOURCES} ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  bcstrategy_composite_factory
  SOURCES bcstrategy_composite_factory.cpp ${UNIT_TEST_DRIVER} user_app_BCStrategy_Dirichlet_Constant.hpp user_app_BCStrategy_Dirichlet_Constant_impl.hpp user_app_BCStrategy_Factory_Physics1.hpp user_app_BCStrategy_Factory_Physics2.hpp
  NUM_MPI_PROCS 1
  )
