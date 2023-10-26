

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_INCLUDE_DIRECTORIES(${PARENT_PACKAGE_SOURCE_DIR}/disc-fe/test/equation_set)
TRIBITS_INCLUDE_DIRECTORIES(${PARENT_PACKAGE_SOURCE_DIR}/disc-fe/test/closure_model)
TRIBITS_INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/test/bcstrategy)

SET(UNIT_TEST_DRIVER
  ${PANZER_UNIT_TEST_MAIN})

IF(NOT Phalanx_KOKKOS_DEVICE_TYPE STREQUAL "CUDA")
  # the structure of this problem combined with repeated reinitialization
  # is causing real problems on CUDA... prevent it!
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    assembly_engine
    SOURCES assembly_engine.cpp ${UNIT_TEST_DRIVER}
    COMM serial mpi
    NUM_MPI_PROCS 4
    )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  simple_bc
  SOURCES simple_bc.cpp ${UNIT_TEST_DRIVER}
  COMM serial mpi
  NUM_MPI_PROCS 2
  )
  
