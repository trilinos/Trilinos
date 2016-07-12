
IF(${PROJECT_NAME}_ENABLE_Epetra)
#  TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    EpetraBatchManagerInterface
#    SOURCES test_01.cpp
#    ARGS PrintItAll
#    COMM serial mpi
#    PASS_REGULAR_EXPRESSION "TEST PASSED"
#    ADD_DIR_TO_NAME
#    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    EpetraSROMSampleGenerator
    SOURCES test_02.cpp
    ARGS PrintItAll
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(
      solEpetraDataCopy
      SOURCE_FILES
#        input_01.xml
        input_02.xml
      SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
      DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
    )
ENDIF()

