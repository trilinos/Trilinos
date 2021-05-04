
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

TRIBITS_ADD_TEST(
  perf_test
  NAME perf_test_a
  ARGS "101551 1"
  COMM mpi
  NUM_MPI_PROCS 12
  RUN_SERIAL
  CATEGORIES PERFORMANCE
  )

TRIBITS_ADD_TEST(
  perf_test
  NAME perf_test_b
  ARGS "101551 2"
  COMM mpi
  NUM_MPI_PROCS 12
  RUN_SERIAL
  CATEGORIES PERFORMANCE
  )

TRIBITS_ADD_TEST(
  perf_test
  NAME perf_test_c
  ARGS "101551 3"
  COMM mpi
  NUM_MPI_PROCS 12
  RUN_SERIAL
  CATEGORIES PERFORMANCE
  )

TRIBITS_ADD_TEST(
  perf_test
  NAME perf_test_d
  ARGS "101551 4"
  COMM mpi
  NUM_MPI_PROCS 12
  RUN_SERIAL
  CATEGORIES PERFORMANCE
  )

TRIBITS_ADD_TEST(
  perf_test
  NAME perf_test_e
  ARGS "101551 6"
  COMM mpi
  NUM_MPI_PROCS 12
  RUN_SERIAL
  CATEGORIES PERFORMANCE
  )

TRIBITS_ADD_TEST(
  perf_test
  NAME perf_test_f
  ARGS "101551 12"
  COMM mpi
  NUM_MPI_PROCS 12
  RUN_SERIAL
  CATEGORIES PERFORMANCE
  )

IF(${PROJECT_NAME}_TEST_CATEGORIES STREQUAL PERFORMANCE)
  ADD_TEST( Adelus_perf_test_MPI_1_M jsrun --smpiargs="-gpu" -n 1 -a 1 -c 10 -g 1 -K 1 -d packed -b packed:10 ./Adelus_perf_test.exe 27882 1 )
  ADD_TEST( Adelus_perf_test_MPI_4_M jsrun --smpiargs="-gpu" -n 4 -a 1 -c 10 -g 1 -K 2 -d packed -b packed:10 ./Adelus_perf_test.exe 55764 2 )
  ADD_TEST( Adelus_perf_test_MPI_9_M jsrun --smpiargs="-gpu" -n 9 -a 1 -c 10 -g 1 -K 2 -d packed -b packed:10 ./Adelus_perf_test.exe 83646 3 )
  ADD_TEST( Adelus_perf_test_MPI_16_M jsrun --smpiargs="-gpu" -n 16 -a 1 -c 10 -g 1 -K 2 -d packed -b packed:10 ./Adelus_perf_test.exe 111528 4 )
  ADD_TEST( Adelus_perf_test_MPI_25_M jsrun --smpiargs="-gpu" -n 25 -a 1 -c 10 -g 1 -K 2 -d packed -b packed:10 ./Adelus_perf_test.exe 139410 5 )

  ADD_TEST( Adelus_perf_test_a_MPI_12_M jsrun --smpiargs="-gpu" -n 12 -a 1 -c 10 -g 1 -K 2 -d packed -b packed:10 ./Adelus_perf_test.exe 101551 1 )
  ADD_TEST( Adelus_perf_test_b_MPI_12_M jsrun --smpiargs="-gpu" -n 12 -a 1 -c 10 -g 1 -K 2 -d packed -b packed:10 ./Adelus_perf_test.exe 101551 2 )
  ADD_TEST( Adelus_perf_test_c_MPI_12_M jsrun --smpiargs="-gpu" -n 12 -a 1 -c 10 -g 1 -K 2 -d packed -b packed:10 ./Adelus_perf_test.exe 101551 3 )
  ADD_TEST( Adelus_perf_test_d_MPI_12_M jsrun --smpiargs="-gpu" -n 12 -a 1 -c 10 -g 1 -K 2 -d packed -b packed:10 ./Adelus_perf_test.exe 101551 4 )
  ADD_TEST( Adelus_perf_test_e_MPI_12_M jsrun --smpiargs="-gpu" -n 12 -a 1 -c 10 -g 1 -K 2 -d packed -b packed:10 ./Adelus_perf_test.exe 101551 6 )
  ADD_TEST( Adelus_perf_test_f_MPI_12_M jsrun --smpiargs="-gpu" -n 12 -a 1 -c 10 -g 1 -K 2 -d packed -b packed:10 ./Adelus_perf_test.exe 101551 12 )

  #ADD_TEST( Adelus_perf_test_MPI_T2_1 jsrun --smpiargs="-gpu" -n 1 -p 1 -c ALL_CPUS -g ALL_GPUS -d packed -b packed:10 ./Adelus_perf_test.exe 27882 1 )
  #ADD_TEST( Adelus_perf_test_MPI_T2_4 jsrun --smpiargs="-gpu" -n 1 -p 4 -c ALL_CPUS -g ALL_GPUS -d packed -b packed:10 ./Adelus_perf_test.exe 55764 2 )
  #ADD_TEST( Adelus_perf_test_MPI_T2_9 jsrun --smpiargs="-gpu" -n 3 -p 9 -c ALL_CPUS -g ALL_GPUS -d packed -b packed:10 ./Adelus_perf_test.exe 83646 3 )
  #ADD_TEST( Adelus_perf_test_MPI_T2_16 jsrun --smpiargs="-gpu" -n 4 -p 16 -c ALL_CPUS -g ALL_GPUS -d packed -b packed:10 ./Adelus_perf_test.exe 111528 4 )
  #ADD_TEST( Adelus_perf_test_MPI_T2_25 jsrun --smpiargs="-gpu" -n 7 -p 25 -c ALL_CPUS -g ALL_GPUS -d packed -b packed:10 ./Adelus_perf_test.exe 139410 5 )

  #ADD_TEST( Adelus_perf_test_a_MPI_T2_12 jsrun --smpiargs="-gpu" -n 3 -p 12 -c ALL_CPUS -g ALL_GPUS -d packed -b packed:10 ./Adelus_perf_test.exe 101551 1 )
  #ADD_TEST( Adelus_perf_test_b_MPI_T2_12 jsrun --smpiargs="-gpu" -n 3 -p 12 -c ALL_CPUS -g ALL_GPUS -d packed -b packed:10 ./Adelus_perf_test.exe 101551 2 )
  #ADD_TEST( Adelus_perf_test_c_MPI_T2_12 jsrun --smpiargs="-gpu" -n 3 -p 12 -c ALL_CPUS -g ALL_GPUS -d packed -b packed:10 ./Adelus_perf_test.exe 101551 3 )
  #ADD_TEST( Adelus_perf_test_d_MPI_T2_12 jsrun --smpiargs="-gpu" -n 3 -p 12 -c ALL_CPUS -g ALL_GPUS -d packed -b packed:10 ./Adelus_perf_test.exe 101551 4 )
  #ADD_TEST( Adelus_perf_test_e_MPI_T2_12 jsrun --smpiargs="-gpu" -n 3 -p 12 -c ALL_CPUS -g ALL_GPUS -d packed -b packed:10 ./Adelus_perf_test.exe 101551 6 )
  #ADD_TEST( Adelus_perf_test_f_MPI_T2_12 jsrun --smpiargs="-gpu" -n 3 -p 12 -c ALL_CPUS -g ALL_GPUS -d packed -b packed:10 ./Adelus_perf_test.exe 101551 12 )
ENDIF()
