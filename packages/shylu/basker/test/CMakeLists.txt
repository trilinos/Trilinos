IF(Kokkos_ENABLE_OpenMP)

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
        amesos2_interface_test
	NOEXEPREFIX
	SOURCES amesos2_interface_test.cpp test_util.hpp
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

ENDIF() ##END Kokkos_ENABLE_OpenMP

TRIBITS_COPY_FILES_TO_BINARY_DIR(${PACKAGE_NAME}_Test_Files
        SOURCE_FILES matrix1.mtx 
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
        DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

                
