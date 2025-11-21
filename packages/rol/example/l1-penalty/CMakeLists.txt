tribits_include_directories(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/vector)
tribits_include_directories(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/function)

ROL_ADD_EXECUTABLE_AND_TEST(
  example_01
  SOURCES example_01.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
)

ROL_COPY_FILES_TO_BINARY_DIR(
  L1PenaltyDataCopy
  SOURCE_FILES
    input.xml
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}" 
)
