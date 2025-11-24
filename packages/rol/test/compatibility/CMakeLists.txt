
ROL_ADD_EXECUTABLE_AND_TEST(
  stacktrace 
  SOURCES test_01.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

add_subdirectory(eigen)

