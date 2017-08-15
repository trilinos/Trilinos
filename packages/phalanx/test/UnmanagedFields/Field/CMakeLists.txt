TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  unmanaged_fields
  SOURCES UnmanagedFieldTest.cpp Field_TestEvaluators.hpp Field_TestEvaluators_Def.hpp
  TESTONLYLIBS phalanx_unit_test_main
  COMM serial mpi
  NUM_MPI_PROCS 1
)

