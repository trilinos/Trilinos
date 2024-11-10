# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Belos
  ${PACKAGE_NAME}_ENABLE_Epetra
  ${PACKAGE_NAME}_ENABLE_Zoltan2Core
  )

#
# Executable
#


#TODO: how to remove that?
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)



TRIBITS_ADD_EXECUTABLE(
  ParameterListInterpreter
  SOURCES ParameterListInterpreter.cpp ../unit_tests/MueLu_TestHelpers.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  CreateOperator
  SOURCES CreateOperator.cpp ../unit_tests/MueLu_TestHelpers.cpp
  COMM serial mpi
  )



#
# Tests
#

# Run the ParameterListInterpreter test with a single thread to avoid non-deterministic behavior
SET(${PACKAGE_NAME}_PARAMETERLISTINTERPRETERTEST_SINGLETHREAD " --kokkos-num-threads=1")


# Test that demangling is enabled, otherwise the test would fail even is they are logically correct.
# NOTE:
# A the moment, the test environment for Trilinos on software is set incorrectly, and
# results in not setting the value correctly.
IF ( ( HAVE_GCC_ABI_DEMANGLE STREQUAL "1") OR ( HAVE_GCC_ABI_DEMANGLE STREQUAL "ON") )

  IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND
      ${PACKAGE_NAME}_ENABLE_Ifpack AND ${PACKAGE_NAME}_ENABLE_Amesos AND
      HAVE_MUELU_INST_DOUBLE_INT_INT)

    string(CONCAT TEST_ARGS "--linAlgebra=Epetra" "${${PACKAGE_NAME}_PARAMETERLISTINTERPRETERTEST_SINGLETHREAD}")

    TRIBITS_ADD_TEST(
      CreateOperator
      NAME "CreateOperatorEpetra"
      ARGS ${TEST_ARGS}
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      COMM mpi serial
      )
    TRIBITS_ADD_TEST(
      CreateOperator
      NAME "CreateOperatorEpetra"
      ARGS ${TEST_ARGS}
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 4
      COMM mpi
      )

    IF(${PACKAGE_NAME}_ENABLE_Zoltan2Core)

      TRIBITS_ADD_TEST(
        ParameterListInterpreter
        NAME "ParameterListInterpreterEpetra"
        ARGS ${TEST_ARGS}
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        COMM mpi serial
      )
      TRIBITS_ADD_TEST(
        ParameterListInterpreter
        NAME "ParameterListInterpreterEpetra"
        ARGS ${TEST_ARGS}
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 4
        COMM mpi
      )

      string(CONCAT TEST_ARGS "--linAlgebra=Epetra --heavytests" "${${PACKAGE_NAME}_PARAMETERLISTINTERPRETERTEST_SINGLETHREAD}")

      # This test only runs on the dashboard, not at checkin.
      TRIBITS_ADD_TEST(
        ParameterListInterpreter
        NAME "ParameterListInterpreterEpetraHeavy"
        ARGS ${TEST_ARGS}
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        COMM mpi serial
        CATEGORIES CONTINUOUS
      )

      # This test only runs on the dashboard, not at checkin.
      TRIBITS_ADD_TEST(
        ParameterListInterpreter
        NAME "ParameterListInterpreterEpetraHeavy"
        ARGS ${TEST_ARGS}
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 4
        COMM mpi
        CATEGORIES CONTINUOUS
      )

    ENDIF()

  ENDIF()

  IF(${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)

    string(CONCAT TEST_ARGS "--linAlgebra=Tpetra" "${${PACKAGE_NAME}_PARAMETERLISTINTERPRETERTEST_SINGLETHREAD}")

    TRIBITS_ADD_TEST(
      CreateOperator
      NAME "CreateOperatorTpetra"
      ARGS ${TEST_ARGS}
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      COMM mpi serial
      )
    TRIBITS_ADD_TEST(
      CreateOperator
      NAME "CreateOperatorTpetra"
      ARGS ${TEST_ARGS}
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 4
      COMM mpi
      )

    IF(${PACKAGE_NAME}_ENABLE_Zoltan2Core)

      TRIBITS_ADD_TEST(
        ParameterListInterpreter
        NAME "ParameterListInterpreterTpetra"
        ARGS ${TEST_ARGS}
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        COMM mpi serial
      )
      TRIBITS_ADD_TEST(
        ParameterListInterpreter
        NAME "ParameterListInterpreterTpetra"
        ARGS ${TEST_ARGS}
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 4
        COMM mpi
      )

      string(CONCAT TEST_ARGS "--linAlgebra=Tpetra --heavytests" "${${PACKAGE_NAME}_PARAMETERLISTINTERPRETERTEST_SINGLETHREAD}")

      # This test only runs on the dashboard, not at checkin.
      TRIBITS_ADD_TEST(
        ParameterListInterpreter
        NAME "ParameterListInterpreterTpetraHeavy"
        ARGS ${TEST_ARGS}
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        COMM mpi serial
        CATEGORIES CONTINUOUS
      )

      # This test only runs on the dashboard, not at checkin.
      TRIBITS_ADD_TEST(
        ParameterListInterpreter
        NAME "ParameterListInterpreterTpetraHeavy"
        ARGS ${TEST_ARGS}
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 4
        COMM mpi
        CATEGORIES CONTINUOUS
      )

    ENDIF()

  ENDIF()

ENDIF()



# Select appropriate set of input and gold files, add make target for rebase.
file (REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/default)
ADD_SUBDIRECTORY(default)

file (REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/kokkos)
ADD_SUBDIRECTORY(kokkos)


TRIBITS_COPY_FILES_TO_BINARY_DIR(rebase_cp
  SOURCE_FILES rebase.sh
  )

ADD_CUSTOM_TARGET(
  rebase
  COMMAND ./rebase.sh "${CMAKE_CURRENT_SOURCE_DIR}"
  )
