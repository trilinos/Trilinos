
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${TEUCHOS_STD_UNIT_TEST_MAIN})

SET(SOURCES
    equation_set.cpp 
    user_app_EquationSetFactory.hpp 
    user_app_EquationSet_Energy.hpp
    user_app_EquationSet_Energy_impl.hpp
    user_app_Convection.hpp
    user_app_Convection_impl.hpp
    ) 

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  equation_set
  SOURCES ${SOURCES} ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
