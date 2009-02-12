
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
  test_tpi_c
  COMM serial
  SOURCES test_main.c test_c_tpi.c test_c_dnax.c test_tpi_unit.c
  DEPLIBS pthread
  DIRECTORY .
  )

PACKAGE_ADD_EXECUTABLE(
  test_tpi_cpp
  COMM serial
  SOURCES test_tpi.cpp
  DEPLIBS pthread
  DIRECTORY .
  )

PACKAGE_ADD_EXECUTABLE(
  test_tpi_sum
  COMM serial mpi
  SOURCES test_mpi_sum.c
  DEPLIBS pthread
  DIRECTORY .
  )

PACKAGE_ADD_TEST(
  test_tpi_c
  NAME test_tpi_c
  ARGS "1"
  COMM serial
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
  ARGS "1 2 4 6 8"
  COMM serial
  DIRECTORY .
  )

PACKAGE_ADD_TEST(
  test_tpi_sum
  NAME test_tpi_sum_np1
  ARGS "1"
  COMM mpi
  NUM_MPI_PROCS 1
  DIRECTORY .
  )

PACKAGE_ADD_TEST(
  test_tpi_sum
  NAME test_tpi_sum_np2
  ARGS "1"
  COMM mpi
  NUM_MPI_PROCS 2
  DIRECTORY .
  )

PACKAGE_ADD_TEST(
  test_tpi_sum
  NAME test_tpi_sum_np4
  ARGS "1"
  COMM mpi
  NUM_MPI_PROCS 4
  DIRECTORY .
  )


