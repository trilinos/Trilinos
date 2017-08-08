TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  unmanaged_mdfields
  SOURCES UnmanagedMDFieldTest.cpp MDField_TestEvaluators.hpp MDField_TestEvaluators_Def.hpp
  TESTONLYLIBS phalanx_unit_test_main
  COMM serial mpi
  NUM_MPI_PROCS 1
)

