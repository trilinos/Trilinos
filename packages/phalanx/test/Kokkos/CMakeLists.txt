TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../Utilities)

# RUN_SERIAL is added since UniqueToken can require a large amount of
# memory on GPUs.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tKokkos
  SOURCES tKokkos.cpp
  TESTONLYLIBS phalanx_unit_test_main phalanx_test_utilities
  NUM_MPI_PROCS 1
  RUN_SERIAL
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

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tKokkosVirtualFunctionOnDevice
  SOURCES tKokkosVirtualFunctionOnDevice.cpp
  TESTONLYLIBS phalanx_unit_test_main phalanx_test_utilities
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tKokkosViewOfViews
  SOURCES tKokkosViewOfViews.cpp
  TESTONLYLIBS phalanx_unit_test_main phalanx_test_utilities
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tKokkosClassOnDevice
  SOURCES tKokkosClassOnDevice.cpp
  TESTONLYLIBS phalanx_unit_test_main phalanx_test_utilities
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tKokkosScalarFunctionsFad
  SOURCES tKokkosScalarFunctionsFad.cpp
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  tKokkosFencing
  SOURCES tKokkosFencing.cpp
  TESTONLYLIBS phalanx_test_utilities
  NUM_MPI_PROCS 1
  )
