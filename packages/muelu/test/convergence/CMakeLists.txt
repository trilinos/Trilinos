# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Belos
  ${PACKAGE_NAME}_ENABLE_Epetra
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
  Convergence
  SOURCES Convergence.cpp ../unit_tests/MueLu_TestHelpers.cpp
  COMM serial mpi
)

#
# Tests
#

IF(${PACKAGE_NAME}_ENABLE_Belos)

  IF(${PACKAGE_NAME}_ENABLE_Epetra
     AND ${PACKAGE_NAME}_ENABLE_EpetraExt
     AND ${PACKAGE_NAME}_ENABLE_Ifpack
     AND ${PACKAGE_NAME}_ENABLE_Amesos
  )

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

  IF(${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)

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
