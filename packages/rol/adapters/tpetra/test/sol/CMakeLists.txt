
IF(${PROJECT_NAME}_ENABLE_Tpetra)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TpetraSimulatedConstraintInterface
    SOURCES test_01.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TpetraSimulatedConstraintInterfaceCVaR
    SOURCES test_02.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TpetraSimulatedConstraintInterfaceDual
    SOURCES test_03.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    TpetraSimulatedConstraintInterfaceDataCopy
    SOURCE_FILES
      input.xml
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )

ENDIF()


