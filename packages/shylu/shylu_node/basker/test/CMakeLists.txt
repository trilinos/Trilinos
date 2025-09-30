#IF(Kokkos_ENABLE_OPENMP)

    TRIBITS_ADD_EXECUTABLE(
        simple_test_check
        NOEXEPREFIX
        SOURCES simple_test_check.cpp
     )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
        singular_matrix_check
        NOEXEPREFIX
        SOURCES singular_matrix_check.cpp
        NUM_MPI_PROCS 1
        FAIL_REGULAR_EXPRESSION "  FAILED  "
     )

    TRIBITS_ADD_EXECUTABLE(
        basker_test
        NOEXEPREFIX
        SOURCES basker_test.cpp
     )

     TRIBITS_ADD_EXECUTABLE(
        basker_test_inc
        NOEXEPREFIX
        SOURCES basker_test_inc.cpp
     )

     TRIBITS_ADD_EXECUTABLE(
        basker_amesos2_interface_test
        NOEXEPREFIX
        SOURCES amesos2_interface_test.cpp test_util.hpp
     )
     
     TRIBITS_ADD_EXECUTABLE(
        basker_amesos2_interface_coverage_test
        NOEXEPREFIX
        SOURCES amesos2_interface_coverage_test.cpp test_util.hpp
     )
         
     ##Single threaded test
     TRIBITS_ADD_ADVANCED_TEST(
         basker_test_1
         TEST_0 EXEC basker_test
             NOEXEPREFIX
              ARGS matrix1.mtx 1
          FAIL_REGULAR_EXPRESSION "FAILED"
          COMM serial mpi
          OVERALL_NUM_MPI_PROCS 1
     )

     TRIBITS_ADD_ADVANCED_TEST(
          basker_amesos_test_1
          TEST_0 EXEC basker_amesos2_interface_test
              NOEXEPREFIX
              ARGS 2 matrix1.mtx
           FAILED_REGULAR_EXPRESSION "FAILED"
           COMM serial mpi
           OVERALL_NUM_MPI_PROCS 1
     )
        
     ##4 threaded test
     TRIBITS_ADD_ADVANCED_TEST(
         basker_test_4
         TEST_0 EXEC basker_test
             NOEXEPREFIX
              ARGS matrix1.mtx 4
          FAIL_REGULAR_EXPRESSION "FAILED"
          COMM serial mpi
          OVERALL_NUM_MPI_PROCS 1
     )
     
     TRIBITS_ADD_ADVANCED_TEST(
         basker_amesos2_test_4
         TEST_0 EXEC basker_amesos2_interface_test
             NOEXEPREFIX
             ARGS 4 matrix1.mtx
          FAIL_REGULAR_EXPRESSION "FAILED"
          COMM serial mpi
          OVERALL_NUM_MPI_PROCS 1
     )

     ##1-4 Coverage Test
     TRIBITS_ADD_ADVANCED_TEST(
          basker_amesos2_coverage_test_14
          TEST_0 EXEC basker_amesos2_interface_coverage_test
               NOEXEPREFIX
          FAIL_REGULAR_EXPRESSION "FAILED"
          COMM serial mpi
          OVERALL_NUM_MPI_PROCS 1
     )

#ENDIF() ##END Kokkos_ENABLE_OPENMP

TRIBITS_COPY_FILES_TO_BINARY_DIR(${PACKAGE_NAME}_Test_Files
        SOURCE_FILES matrix1.mtx coverage.mtx 
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
        DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)


