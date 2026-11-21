TRIBITS_ADD_EXECUTABLE(
  test_generate_len_name32
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES test_generate_len_name32.c
  COMM serial mpi
  )

IF (${CMAKE_PROJECT_NAME}_ENABLE_SEACASExodus)

TRIBITS_ADD_ADVANCED_TEST(
 exodiff_test_len_name_32
 TEST_0 EXEC ${CMAKE_CURRENT_BINARY_DIR}/test_generate_len_name32
        NOEXEPREFIX NOEXESUFFIX
 TEST_1 CMND ${CMAKE_CURRENT_BINARY_DIR}/../exodiff
        ARGS -f ${CMAKE_CURRENT_SOURCE_DIR}/len-name-32.exodiff
		         ${CMAKE_CURRENT_BINARY_DIR}/len-name-32-a.exo
		         ${CMAKE_CURRENT_BINARY_DIR}/len-name-32-b.exo
		    OUTPUT_FILE exodiff-len-name-32.stdout
		    PASS_REGULAR_EXPRESSION "exodiff: Files are the same"
 TEST_2 CMND ${CMAKE_CURRENT_BINARY_DIR}/../exodiff
        ARGS -f ${CMAKE_CURRENT_SOURCE_DIR}/len-name-32.exodiff
		         ${CMAKE_CURRENT_BINARY_DIR}/len-name-32-a.exo
		         ${CMAKE_CURRENT_BINARY_DIR}/len-name-32-b.exo
		         ${CMAKE_CURRENT_BINARY_DIR}/len-name-32-diff.exo
		    OUTPUT_FILE exodiff-len-name-32-diff.stdout
		    FAIL_REGULAR_EXPRESSION "is too long"
 COMM mpi serial
 OVERALL_NUM_MPI_PROCS 1
 FINAL_PASS_REGULAR_EXPRESSION
 XHOSTTYPE Windows
)

endif()
