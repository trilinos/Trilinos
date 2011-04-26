
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageLibraryMacros)

####################

SET(HEADERS "")
SET(SOURCES "")
 
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
 
SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )
 
INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
 
APPEND_SET(HEADERS
  BoxPartition.h
  CGSolver.h
  tpi_vector.h
  dcrs_matrix.h
  )
 
####################


PACKAGE_ADD_EXECUTABLE(
  test_tpi_hhpccg
  COMM serial mpi
  SOURCES main.c CGSolver.c BoxPartitionIB.c tpi_vector.c dcrs_matrix.c
  DEPLIBS pthread m
  DIRECTORY .
  )

PACKAGE_ADD_TEST(
  test_tpi_hhpccg
  NAME test_tpi_hhpccg_serial_1
  COMM serial
  DIRECTORY .
  )

PACKAGE_ADD_TEST(
  test_tpi_hhpccg
  NAME test_tpi_hhpccg_serial_2
  COMM serial
  ARGS "threads=2"
  DIRECTORY .
  )

PACKAGE_ADD_TEST(
  test_tpi_hhpccg
  NAME test_tpi_hhpccg_serial_4
  COMM serial
  ARGS "threads=4"
  DIRECTORY .
  XHOSTTYPE AIX
  )

PACKAGE_ADD_TEST(
  test_tpi_hhpccg
  NAME test_tpi_hhpccg_mpi_1
  COMM mpi
  NUM_MPI_PROCS 1
  DIRECTORY .
  )

PACKAGE_ADD_TEST(
  test_tpi_hhpccg
  NAME test_tpi_hhpccg_mpi_2
  COMM mpi
  NUM_MPI_PROCS 2
  DIRECTORY .
  )

PACKAGE_ADD_TEST(
  test_tpi_hhpccg
  NAME test_tpi_hhpccg_mpi_4
  COMM mpi
  NUM_MPI_PROCS 4
  DIRECTORY .
  )


