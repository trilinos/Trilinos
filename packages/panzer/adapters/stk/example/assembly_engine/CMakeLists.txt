INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${TRILINOS_HOME_DIR}/Panzer/panzer/test/equation_set)

INCLUDE_DIRECTORIES(${TRILINOS_HOME_DIR}/Panzer/panzer/test/model_factory)

INCLUDE_DIRECTORIES(${TRILINOS_HOME_DIR}/Panzer/panzer/test/bcstrategy)
  
SET(ASSEMBLY_EXAMPLE_SOURCES
  assembly_driver.cpp
  write_solution_data.cpp
  )
  
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  assembly_example
  SOURCES ${ASSEMBLY_EXAMPLE_SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
