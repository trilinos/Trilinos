
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
  test_tpi_c
  COMM serial
  SOURCES test_main.c test_c_dnax.c test_tpi_unit.c
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
  COMM serial
  SOURCES test_mpi_sum.c
  DEPLIBS pthread
  DIRECTORY .
  )

PACKAGE_ADD_EXECUTABLE(
  test_tpi_sum_mpi
  COMM mpi
  SOURCES test_mpi_sum.c
  DEPLIBS pthread
  DIRECTORY .
  )

PACKAGE_ADD_TEST(
  test_tpi_c
  NAME test_tpi_c
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
  COMM serial
  DIRECTORY .
  )

PACKAGE_ADD_TEST(
  test_tpi_sum_mpi
  NAME test_tpi_sum_np1
  COMM mpi
  NUM_MPI_PROCS 1
  DIRECTORY .
  )

PACKAGE_ADD_TEST(
  test_tpi_sum_mpi
  NAME test_tpi_sum_np2
  COMM mpi
  NUM_MPI_PROCS 2
  DIRECTORY .
  )

PACKAGE_ADD_TEST(
  test_tpi_sum_mpi
  NAME test_tpi_sum_np4
  COMM mpi
  NUM_MPI_PROCS 4
  DIRECTORY .
  )


