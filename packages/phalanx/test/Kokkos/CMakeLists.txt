INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../Utilities)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tKokkos
  SOURCES tKokkos.cpp
  TESTONLYLIBS phalanx_unit_test_main phalanx_test_utilities
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tKokkosPerf
  SOURCES tKokkosPerf.cpp
  TESTONLYLIBS phalanx_unit_test_main phalanx_test_utilities
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tKokkosNestedLambda
  SOURCES tKokkosNestedLambda.cpp
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tKokkosAllocationTracker
  SOURCES tKokkosAllocationTracker.cpp
  TESTONLYLIBS phalanx_unit_test_main phalanx_test_utilities
  NUM_MPI_PROCS 1
  )
