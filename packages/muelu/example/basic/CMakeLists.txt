INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Executables
#

#
#Simple-LongInt
#
IF(${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)
  IF(Teuchos_ENABLE_LONG_LONG_INT)

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        Simple-LongInt
        SOURCES Simple-LongInt.cpp
        DEPLIBS muelu xpetra
        NUM_MPI_PROCS 4
        COMM serial mpi
        )

  ENDIF()
ENDIF()

#
#Simple-Complex
#
IF(${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)

  IF(HAVE_MUELU_COMPLEX)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Simple-Complex
      SOURCES Simple-Complex.cpp
      DEPLIBS muelu xpetra
      NUM_MPI_PROCS 1
      COMM serial mpi
      )
  ENDIF()

ENDIF()

#TODO: NUM_MPI_PROCS->4
