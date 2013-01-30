INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER
  ${TEUCHOS_STD_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  integration_rule
  SOURCES integration_rule.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  integration_values
  SOURCES integration_values.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  dimension
  SOURCES dimension.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  version
  SOURCES version.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  basis
  SOURCES basis.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  basis_values
  SOURCES basis_values.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  point_values 
  SOURCES point_values.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  boundary_condition
  SOURCES boundary_condition.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  material_model_entry
  SOURCES material_model_entry.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  stlmap_utilities
  SOURCES stlmap_utilities.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  shards_utilities
  SOURCES shards_utilities.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  string_utilities
  SOURCES string_utilities.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  evaluators
  SOURCES evaluators.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  element_block_to_physics_block_map
  SOURCES element_block_to_physics_block_map.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  zero_sensitivities
  SOURCES zero_sensitivities.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  output_stream
  SOURCES output_stream.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  global_data
  SOURCES global_data.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  parameter_library
  SOURCES parameter_library.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  cell_topology_info
  SOURCES cell_topology_info.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  parameter_list_acceptance_test
  SOURCES parameter_list_acceptance_test.cpp ${UNIT_TEST_DRIVER}
  NUM_MPI_PROCS 1
  )

IF(Panzer_BUILD_PAPI_SUPPORT)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    papi_raw
    SOURCES papi_raw.cpp ${UNIT_TEST_DRIVER}
    NUM_MPI_PROCS 1
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    papi_decorator
    SOURCES papi_decorator.cpp ${UNIT_TEST_DRIVER}
    NUM_MPI_PROCS 1
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    papi_improved
    SOURCES papi_improved.cpp ${UNIT_TEST_DRIVER}
    NUM_MPI_PROCS 1
    )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(core_tests
  SOURCE_FILES parameter_list_acceptance_test.xml parameter_list_acceptance_test2.xml
  EXEDEPS parameter_list_acceptance_test
  )
