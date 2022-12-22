TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/sol)
TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/mpi/src/vector)

IF(${PROJECT_NAME}_ENABLE_Tpetra)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TpetraMultiVectorInterface
    SOURCES test_01.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )
ENDIF()

IF(${PROJECT_NAME}_ENABLE_Tpetra)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SimulatedVectorTpetraBatchManagerInterface
    SOURCES test_02.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )
ENDIF()

IF(${PROJECT_NAME}_ENABLE_Tpetra)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ScaledTpetraMultiVector
    SOURCES test_03.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )
ENDIF()

IF(${PROJECT_NAME}_ENABLE_Tpetra)
#  TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    PinTTpetraMultiVector
#    SOURCES test_04.cpp
#    ARGS PrintItAll
#    COMM mpi
#    NUM_MPI_PROCS 3
#    PASS_REGULAR_EXPRESSION "TEST PASSED"
#    ADD_DIR_TO_NAME
#    )
ENDIF()
