INCLUDE(TribitsAddExecutableAndTest)

IF(${PROJECT_NAME}_ENABLE_Belos)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    BelosInterface
    SOURCES test_01.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    ) 

    TRIBITS_COPY_FILES_TO_BINARY_DIR(
    BelosInterfaceDataCopy
    SOURCE_FILES
     parameters.xml
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )
ENDIF()
