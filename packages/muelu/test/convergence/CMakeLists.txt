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

INCLUDE(TribitsAddExecutable)

#TODO: how to remove that?
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

IF (${PACKAGE_NAME}_ENABLE_Galeri)

  TRIBITS_ADD_EXECUTABLE(
    Convergence
    SOURCES Convergence.cpp ../unit_tests/MueLu_TestHelpers.cpp
    COMM serial mpi
    )

ENDIF()

#
# Tests
#

IF (${PACKAGE_NAME}_ENABLE_Galeri AND ${PACKAGE_NAME}_ENABLE_Belos)

  IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND
      ${PACKAGE_NAME}_ENABLE_Ifpack AND ${PACKAGE_NAME}_ENABLE_Amesos)

    TRIBITS_ADD_TEST(
      Convergence
      NAME "ConvergenceEpetra"
      ARGS "--linAlgebra=Epetra"
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      COMM mpi serial
      )
    TRIBITS_ADD_TEST(
      Convergence
      NAME "ConvergenceEpetra"
      ARGS "--linAlgebra=Epetra"
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 4
      COMM mpi
      )

  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)

    TRIBITS_ADD_TEST(
      Convergence
      NAME "ConvergenceTpetra"
      ARGS "--linAlgebra=Tpetra"
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      COMM mpi serial
      )
    TRIBITS_ADD_TEST(
      Convergence
      NAME "ConvergenceTpetra"
      ARGS "--linAlgebra=Tpetra"
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 4
      COMM mpi
      )

  ENDIF()

ENDIF()

ADD_SUBDIRECTORY(Convergence)
