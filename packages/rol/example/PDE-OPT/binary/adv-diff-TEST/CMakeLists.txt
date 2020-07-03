INCLUDE(TribitsAddExecutableAndTest)

IF(${PROJECT_NAME}_ENABLE_Intrepid AND 
   ${PROJECT_NAME}_ENABLE_Ifpack2  AND
   ${PROJECT_NAME}_ENABLE_MueLu    AND
   ${PROJECT_NAME}_ENABLE_Amesos2  AND
   ${PROJECT_NAME}_ENABLE_pebbl    AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  # Need ROL_TpetraMultiVector.hpp
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)

  # Need ROL_TpetraTeuchosBatchManager.hpp
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/sol)

  # Need ROL_Teuchos_LinearOperator.hpp
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/function)

  # Need ROL_TeuchosVector.hpp
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/vector)

  # Need Teuchos/PEBBL interface
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/pebbl)

  # Need ROL/PEBBL interface
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/pebbl/src/algorithm)
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/pebbl/src/function)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    RelaxedBinary
    SOURCES example_01.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
  )

#  TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    SumUpRounding
#    SOURCES example_02.cpp
#    ARGS PrintItAll
#    NUM_MPI_PROCS 4
#    PASS_REGULAR_EXPRESSION "TEST PASSED"
#    ADD_DIR_TO_NAME
#  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    BinaryAdvTESTDiffDataCopy
    SOURCE_FILES
      input_ex01.xml
      sensor.txt
      run.sh
      plotPEBBL.m
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )


ENDIF()
