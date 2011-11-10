INCLUDE(PackageAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Executables
#

IF(( ${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack  AND ${PACKAGE_NAME}_ENABLE_Amesos) OR 
    (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2))
  
  #Simple
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Simple
    SOURCES Simple.cpp
    DEPLIBS muelu muelu-xpetra
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

  #Simple-LongInt
  IF(Teuchos_ENABLE_LONG_LONG_INT)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      Simple-LongInt
      SOURCES Simple-LongInt.cpp
      DEPLIBS muelu muelu-xpetra
      NUM_MPI_PROCS 4
      COMM serial mpi
      )
  ENDIF()

  #Simple-Templated
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Simple-Templated
    SOURCES Simple-Templated.cpp
    DEPLIBS muelu muelu-xpetra
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

  #Simple-Belos
  IF(${PACKAGE_NAME}_ENABLE_Belos)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      Simple-Belos
      SOURCES Simple-Belos.cpp
      NUM_MPI_PROCS 4
      COMM serial mpi
      DEPLIBS muelu muelu-xpetra muelu-adapters
      )
  ENDIF()

  #Simple-Tpetra
  IF(${PACKAGE_NAME}_ENABLE_Tpetra)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      Simple-Tpetra
      SOURCES Simple-Tpetra.cpp
      NUM_MPI_PROCS 4
      COMM serial mpi
      DEPLIBS muelu muelu-xpetra muelu-adapters
      )
  ENDIF()

  #Simple-Epetra
  IF(${PACKAGE_NAME}_ENABLE_Epetra)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      Simple-Epetra
      SOURCES Simple-Epetra.cpp
      NUM_MPI_PROCS 4
      COMM serial mpi
      DEPLIBS muelu muelu-xpetra muelu-adapters
      )
  ENDIF()

ENDIF()
