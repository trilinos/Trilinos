IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_TPETRA)
  # Add current binary dir to include directories
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

  # Create the tests.
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_MultiVectorOpsTests
    SOURCES
        MultiVectorOpsTests.cpp
        ${TEUCHOS_STD_UNIT_TEST_MAIN}
    ARGS ""
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_VectorOpsTests
    SOURCES
        VectorOpsTests.cpp
        ${TEUCHOS_STD_UNIT_TEST_MAIN}
    ARGS ""
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

ENDIF()