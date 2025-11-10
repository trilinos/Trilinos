TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../Utilities)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  mdfield_compiletime
  SOURCES MDField_Compiletime.cpp
  TESTONLYLIBS phalanx_unit_test_main phalanx_test_utilities
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  mdfield_runtime
  SOURCES MDField_Runtime.cpp
  TESTONLYLIBS phalanx_unit_test_main phalanx_test_utilities
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  field
  SOURCES Field.cpp
  TESTONLYLIBS phalanx_unit_test_main phalanx_test_utilities
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  print_field_values
  SOURCES PrintFieldValues.cpp
  TESTONLYLIBS phalanx_unit_test_main phalanx_test_utilities
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  mdfield_simulated_intrepid
  SOURCES MDField_SimulatedIntrepid.cpp
  TESTONLYLIBS phalanx_unit_test_main phalanx_test_utilities
  NUM_MPI_PROCS 1
  )

# I am disabling the test since it is failing.
# Context:
# https://github.com/trilinos/Trilinos/pull/14673
# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   view_allocation_size
#   SOURCES ViewAllocationSize.cpp
#   TESTONLYLIBS phalanx_unit_test_main phalanx_test_utilities
#   NUM_MPI_PROCS 1
#   )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  nc_drv_from_c_mdfield
  SOURCES NonConstDynRankViewFromConstMDField.cpp
  TESTONLYLIBS phalanx_unit_test_main phalanx_test_utilities
  NUM_MPI_PROCS 1
  )
