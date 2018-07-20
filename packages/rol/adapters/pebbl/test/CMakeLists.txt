
IF(ROL_ENABLE_pebbl)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SimpleBinary
    SOURCES test_01.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    pebblConstraint
    SOURCES test_02.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SimpleBinaryEnumeration
    SOURCES test_03.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SimpleBinaryTransformed
    SOURCES test_04.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    pebblDriver
    SOURCES test_05.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    SimpleBinaryDataCopy
    SOURCE_FILES
      input.xml
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )

ENDIF()


