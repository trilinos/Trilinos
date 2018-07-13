IF(${PROJECT_NAME}_ENABLE_Intrepid AND 
   ${PROJECT_NAME}_ENABLE_Epetra   AND
   ${PROJECT_NAME}_ENABLE_Amesos   AND
   ${PROJECT_NAME}_ENABLE_Amesos2  AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/sol)

  # Need ROL_TpetraMultiVector.hpp
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/sol)

  TRIBITS_ADD_EXECUTABLE(
    example_01
    SOURCES example_01.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_ADD_EXECUTABLE(
    example_02
    SOURCES example_02.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_ADD_TEST(
    example_01
    ARGS PrintItAll
    NUM_MPI_PROCS 1
    NUM_TOTAL_CORES_USED 1
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    CATEGORIES HEAVY  # Currently takes 18m to run (see #442)!
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    OldElasticitySIMPOperatorsCopy
    SOURCE_FILES
      input.xml stochastic.xml plotdensity.m plotstate.m plotnodaldensity.m
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )

ENDIF()
