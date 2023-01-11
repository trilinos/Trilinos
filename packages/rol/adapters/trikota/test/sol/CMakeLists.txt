
IF(${PROJECT_NAME}_ENABLE_TriKota)

  TRIBITS_INCLUDE_DIRECTORIES(${ROL_SOURCE_DIR}/adapters/teuchos/src/sol)

  #TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #  SparseGridBurgersControl
  #  SOURCES test_01.cpp
  #  ARGS PrintItAll
  #  COMM serial mpi  
  #  PASS_REGULAR_EXPRESSION "TEST PASSED"
  #  ADD_DIR_TO_NAME
  #  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SparseGridQuadratureTest
    SOURCES test_02.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    SparseGridBurgersControlDataCopy
    SOURCE_FILES
      input.xml
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )

ENDIF()


