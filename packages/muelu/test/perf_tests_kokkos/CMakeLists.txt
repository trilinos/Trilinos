#
# Executable
#

#TODO: how to remove that?
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

IF((${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)
   OR (${PACKAGE_NAME}_ENABLE_Epetra
       AND ${PACKAGE_NAME}_ENABLE_Ifpack
       AND ${PACKAGE_NAME}_ENABLE_Amesos)
)

  TRIBITS_ADD_EXECUTABLE(
    Redirection
    SOURCES Redirection.cpp
    COMM serial mpi
  )

ENDIF()

TRIBITS_ADD_EXECUTABLE(
  KokkosKernels
  SOURCES KokkosKernels.cpp
  COMM serial mpi
)

#
# Tests
#

IF(${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)

  TRIBITS_ADD_TEST(
    Redirection
    NAME "Redirection_Tpetra"
    ARGS "--linAlgebra=Tpetra"
    NUM_MPI_PROCS 4
    COMM mpi # HAVE_MPI required
  )

ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Epetra
   AND ${PACKAGE_NAME}_ENABLE_Ifpack
   AND ${PACKAGE_NAME}_ENABLE_Amesos
)

  TRIBITS_ADD_TEST(
    Redirection
    NAME "Redirection_Epetra"
    ARGS "--linAlgebra=Epetra"
    NUM_MPI_PROCS 4
    COMM mpi # HAVE_MPI required
  )

ENDIF()
