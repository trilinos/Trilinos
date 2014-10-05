
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsLibraryMacros)

####################

SET(HEADERS "")
SET(SOURCES "")
 
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
 
SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )
 
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
 
APPEND_SET(HEADERS
  BoxPartition.h
  CGSolver.h
  tpi_vector.h
  )
 
####################


TRIBITS_ADD_EXECUTABLE(
  test_tpi_hpccg
  COMM serial mpi
  SOURCES main.c CGSolver.c BoxPartition.c tpi_vector.c
  DEPLIBS pthread m
  DIRECTORY .
  )

TRIBITS_ADD_TEST(
  test_tpi_hpccg
  NAME test_tpi_hpccg_serial_1
  COMM serial
  DIRECTORY .
  )

TRIBITS_ADD_TEST(
  test_tpi_hpccg
  NAME test_tpi_hpccg_serial_2
  COMM serial
  ARGS "threads=2"
  DIRECTORY .
  XHOSTTYPE AIX
  )

TRIBITS_ADD_TEST(
  test_tpi_hpccg
  NAME test_tpi_hpccg_serial_4
  COMM serial
  ARGS "threads=4"
  DIRECTORY .
  )

TRIBITS_ADD_TEST(
  test_tpi_hpccg
  NAME test_tpi_hpccg_mpi_1
  COMM mpi
  NUM_MPI_PROCS 1
  DIRECTORY .
  )

TRIBITS_ADD_TEST(
  test_tpi_hpccg
  NAME test_tpi_hpccg_mpi_2
  COMM mpi
  NUM_MPI_PROCS 2
  DIRECTORY .
  )

TRIBITS_ADD_TEST(
  test_tpi_hpccg
  NAME test_tpi_hpccg_mpi_4
  COMM mpi
  NUM_MPI_PROCS 4
  DIRECTORY .
  )


