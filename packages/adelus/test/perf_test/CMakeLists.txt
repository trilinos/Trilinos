
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  perf_test
  SOURCES cxx_main.cpp
  ARGS "27882 1"
  COMM mpi
  NUM_MPI_PROCS 1
  RUN_SERIAL
  CATEGORIES PERFORMANCE
  )

TRIBITS_ADD_TEST(
  perf_test
  ARGS "55764 2"
  COMM mpi
  NUM_MPI_PROCS 4
  RUN_SERIAL
  CATEGORIES PERFORMANCE
  )

TRIBITS_ADD_TEST(
  perf_test
  ARGS "83646 3"
  COMM mpi
  NUM_MPI_PROCS 9
  RUN_SERIAL
  CATEGORIES PERFORMANCE
  )

TRIBITS_ADD_TEST(
  perf_test
  ARGS "111528 4"
  COMM mpi
  NUM_MPI_PROCS 16
  RUN_SERIAL
  CATEGORIES PERFORMANCE
  )

TRIBITS_ADD_TEST(
  perf_test
  ARGS "139410 5"
  COMM mpi
  NUM_MPI_PROCS 25
  RUN_SERIAL
  CATEGORIES PERFORMANCE
  )

IF(${PROJECT_NAME}_TEST_CATEGORIES STREQUAL PERFORMANCE)
  ADD_TEST( Adelus_perf_test_MPI_1_t jsrun --smpiargs="-gpu" -n 1 -a 1 -c 10 -g 1 -K 1 -d packed -b packed:10 ./Adelus_perf_test.exe 27882 1 )

  ADD_TEST( Adelus_perf_test_MPI_4_t jsrun --smpiargs="-gpu" -n 4 -a 1 -c 10 -g 1 -K 2 -d packed -b packed:10 ./Adelus_perf_test.exe 55764 2 )

  ADD_TEST( Adelus_perf_test_MPI_9_t jsrun --smpiargs="-gpu" -n 9 -a 1 -c 10 -g 1 -K 2 -d packed -b packed:10 ./Adelus_perf_test.exe 83646 3 )

  ADD_TEST( Adelus_perf_test_MPI_16_t jsrun --smpiargs="-gpu" -n 16 -a 1 -c 10 -g 1 -K 2 -d packed -b packed:10 ./Adelus_perf_test.exe 111528 4 )

  ADD_TEST( Adelus_perf_test_MPI_25_t jsrun --smpiargs="-gpu" -n 25 -a 1 -c 10 -g 1 -K 2 -d packed -b packed:10 ./Adelus_perf_test.exe 139410 5 )

ENDIF()
