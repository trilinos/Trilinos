IF(Kokkos_ENABLE_SERIAL)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    test_serial
    NOEXEPREFIX
    SOURCES main_serial.cpp 
    FAIL_REGULAR_EXPRESSION "FAILED"
    COMM serial mpi
    NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    test_graphviz
    NOEXEPREFIX
    SOURCES main_graphviz.cpp 
    FAIL_REGULAR_EXPRESSION "FAILED"
    COMM serial mpi
    NUM_MPI_PROCS 1
)
ENDIF()

IF(Kokkos_ENABLE_PTHREAD)
  TRIBITS_ADD_EXECUTABLE(
    test_pthread
    NOEXEPREFIX
    SOURCES main_pthread.cpp 
  )
  INCLUDE(ProcessorCount)
  ProcessorCount(nthreads)
  SET(i 1)
  MATH(EXPR nthreads "${nthreads}+1")
  WHILE(i LESS nthreads)
    TRIBITS_ADD_ADVANCED_TEST(
      test_pthread_${i}
      TEST_0 EXEC test_pthread
             NOEXEPREFIX        
             ARGS --nthreads=${i} 
      FAIL_REGULAR_EXPRESSION "FAILED"
      COMM serial mpi
      OVERALL_NUM_MPI_PROCS 1
    )
    MATH(EXPR i "${i}*2")
  ENDWHILE()
  TRIBITS_ADD_EXECUTABLE(
    test_debug_pthread
    NOEXEPREFIX
    SOURCES main_debug_pthread.cpp 
  )
ENDIF()

IF(Kokkos_ENABLE_QTHREAD)
  TRIBITS_ADD_EXECUTABLE(
    test_qthread
    NOEXEPREFIX
    SOURCES main_qthread.cpp 
  )
  INCLUDE(ProcessorCount)
  ProcessorCount(nthreads)
  SET(i 1)
  MATH(EXPR nthreads "${nthreads}+1")
  WHILE(i LESS nthreads)
    TRIBITS_ADD_ADVANCED_TEST(
      test_qthread_${i}
      TEST_0 EXEC test_qthread
             NOEXEPREFIX        
             ARGS --nthreads=${i} --stack-size=8192
      FAIL_REGULAR_EXPRESSION "FAILED"
      COMM serial mpi
      OVERALL_NUM_MPI_PROCS 1
    )
    MATH(EXPR i "${i}*2")
  ENDWHILE()
  TRIBITS_ADD_EXECUTABLE(
    test_debug_qthread
    NOEXEPREFIX
    SOURCES main_debug_qthread.cpp 
  )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(${PACKAGE_NAME}_Test_Files
        SOURCE_FILES mm_crs_input.mtx mm_crs_chol.mtx
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
        DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
