
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
  test_tpi_unit
  COMM serial mpi
  SOURCES test_tpi_unit.c
  DIRECTORY .
  )

PACKAGE_ADD_EXECUTABLE(
  test_c_dnax
  COMM serial
  SOURCES test_c_dnax.c
  DIRECTORY .
  )

PACKAGE_ADD_EXECUTABLE(
  test_tpi_cpp
  COMM serial
  SOURCES test_tpi.cpp
  DIRECTORY .
  )

PACKAGE_ADD_EXECUTABLE(
  test_tpi_sum
  COMM serial mpi
  SOURCES test_mpi_sum.c
  DIRECTORY .
  )

PACKAGE_ADD_TEST(
  test_tpi_unit
  NAME test_tpi_unit_serial
  COMM serial
  DIRECTORY .
  XHOSTTYPE AIX
  )

PACKAGE_ADD_TEST(
  test_tpi_unit
  NAME test_tpi_unit_mpi
  COMM mpi
  NUM_MPI_PROCS 1
  DIRECTORY .
  )

PACKAGE_ADD_TEST(
  test_tpi_cpp
  NAME test_tpi_cpp
  COMM serial
  DIRECTORY .
  )

PACKAGE_ADD_TEST(
  test_tpi_sum
  NAME test_tpi_sum_serial
  COMM serial
  DIRECTORY .
  XHOSTTYPE AIX
  )

PACKAGE_ADD_TEST(
  test_tpi_sum
  NAME test_tpi_sum_np1
  COMM mpi
  NUM_MPI_PROCS 1
  DIRECTORY .
  )

PACKAGE_ADD_TEST(
  test_tpi_sum
  NAME test_tpi_sum_np2
  COMM mpi
  NUM_MPI_PROCS 2
  DIRECTORY .
  )

PACKAGE_ADD_TEST(
  test_tpi_sum
  NAME test_tpi_sum_np4
  COMM mpi
  NUM_MPI_PROCS 4
  DIRECTORY .
  )


