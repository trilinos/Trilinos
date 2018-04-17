# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Belos
  ${PACKAGE_NAME}_ENABLE_Epetra
  ${PACKAGE_NAME}_ENABLE_Tpetra
  )

#
# Executable
#


#TODO: how to remove that?
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

IF (${PACKAGE_NAME}_ENABLE_Galeri)

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


ENDIF()

#
# Tests
#

IF (${PACKAGE_NAME}_ENABLE_Galeri)

  # Test that demangling is enabled, otherwise the test would fail even is they are logically correct.
  # NOTE:
  # A the moment, the test environment for Trilinos on software is set incorrectly, and
  # results in not setting the value correctly.
  IF ( ( HAVE_GCC_ABI_DEMANGLE STREQUAL "1") OR ( HAVE_GCC_ABI_DEMANGLE STREQUAL "ON") )

    IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND
        ${PACKAGE_NAME}_ENABLE_Ifpack AND ${PACKAGE_NAME}_ENABLE_Amesos AND
        HAVE_MUELU_INST_DOUBLE_INT_INT)

      TRIBITS_ADD_TEST(
        ParameterListInterpreter
        NAME "ParameterListInterpreterEpetra"
        ARGS "--linAlgebra=Epetra"
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        COMM mpi serial
        )
      TRIBITS_ADD_TEST(
        ParameterListInterpreter
        NAME "ParameterListInterpreterEpetra"
        ARGS "--linAlgebra=Epetra"
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 4
        COMM mpi
        )

      TRIBITS_ADD_TEST(
        CreateOperator
        NAME "CreateOperatorEpetra"
        ARGS "--linAlgebra=Epetra"
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        COMM mpi serial
        )
      TRIBITS_ADD_TEST(
        CreateOperator
        NAME "CreateOperatorEpetra"
        ARGS "--linAlgebra=Epetra"
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 4
        COMM mpi
        )

      # This test only runs on the dashboard, not at checkin.
      TRIBITS_ADD_TEST(
        ParameterListInterpreter
        NAME "ParameterListInterpreterEpetraHeavy"
        ARGS "--linAlgebra=Epetra --heavytests"
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        COMM mpi serial
        CATEGORIES CONTINUOUS
        )

      # This test only runs on the dashboard, not at checkin.
      TRIBITS_ADD_TEST(
        ParameterListInterpreter
        NAME "ParameterListInterpreterEpetraHeavy"
        ARGS "--linAlgebra=Epetra --heavytests"
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 4
        COMM mpi
        CATEGORIES CONTINUOUS
        )

    ENDIF()

    IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)

      TRIBITS_ADD_TEST(
        ParameterListInterpreter
        NAME "ParameterListInterpreterTpetra"
        ARGS "--linAlgebra=Tpetra"
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        COMM mpi serial
        )
      TRIBITS_ADD_TEST(
        ParameterListInterpreter
        NAME "ParameterListInterpreterTpetra"
        ARGS "--linAlgebra=Tpetra"
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 4
        COMM mpi
        )

      TRIBITS_ADD_TEST(
        CreateOperator
        NAME "CreateOperatorTpetra"
        ARGS "--linAlgebra=Tpetra"
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        COMM mpi serial
        )
      TRIBITS_ADD_TEST(
        CreateOperator
        NAME "CreateOperatorTpetra"
        ARGS "--linAlgebra=Tpetra"
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 4
        COMM mpi
        )

      # This test only runs on the dashboard, not at checkin.
      TRIBITS_ADD_TEST(
        ParameterListInterpreter
        NAME "ParameterListInterpreterTpetraHeavy"
        ARGS "--linAlgebra=Tpetra --heavytests"
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        COMM mpi serial
        CATEGORIES CONTINUOUS
        )

      # This test only runs on the dashboard, not at checkin.
      TRIBITS_ADD_TEST(
        ParameterListInterpreter
        NAME "ParameterListInterpreterTpetraHeavy"
        ARGS "--linAlgebra=Tpetra --heavytests"
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 4
        COMM mpi
        CATEGORIES CONTINUOUS
        )

    ENDIF()

  ENDIF()

ENDIF()


FILE(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/EasyParameterListInterpreter)
FILE(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/EasyParameterListInterpreter-heavy)
FILE(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/FactoryParameterListInterpreter)
FILE(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/FactoryParameterListInterpreter-heavy)
FILE(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/MLParameterListInterpreter)
FILE(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/MLParameterListInterpreter2)
FILE(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/Output)
IF (NOT ${PACKAGE_NAME}_INST_COMPLEX_INT_INT AND NOT ${PACKAGE_NAME}_INST_COMPLEX_INT_LONG_LONG)
  ADD_SUBDIRECTORY(EasyParameterListInterpreter)
  ADD_SUBDIRECTORY(EasyParameterListInterpreter-heavy)
  ADD_SUBDIRECTORY(FactoryParameterListInterpreter)
  ADD_SUBDIRECTORY(FactoryParameterListInterpreter-heavy)
  ADD_SUBDIRECTORY(MLParameterListInterpreter)
  ADD_SUBDIRECTORY(MLParameterListInterpreter2)
  ADD_SUBDIRECTORY(Output)
ELSE()
  ADD_SUBDIRECTORY(complex)
ENDIF()
