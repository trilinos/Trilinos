
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../Utilities/evaluators)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  performance_allocators
  SOURCES Performance_Allocators.cpp
  DEPLIBS phalanx_test_utilities
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Test passed!"
  )