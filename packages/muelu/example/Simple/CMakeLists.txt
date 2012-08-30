INCLUDE(TribitsAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Executables
#

IF((${PACKAGE_NAME}_ENABLE_Tpetra  AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2) OR 
    (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Ifpack  AND ${PACKAGE_NAME}_ENABLE_Amesos))
  
  #Simple
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Simple
    SOURCES Simple.cpp
    DEPLIBS muelu muelu-xpetra
    NUM_MPI_PROCS 4
    COMM serial mpi
    )
  
  #Simple-LongInt
  IF(Teuchos_ENABLE_LONG_LONG_INT)

    IF(NOT ${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION) #FIXME
      #FIXME no Tpetra instantiation for long long int
      
      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        Simple-LongInt
        SOURCES Simple-LongInt.cpp
        DEPLIBS muelu muelu-xpetra
        NUM_MPI_PROCS 4
        COMM serial mpi
        )
    ENDIF()
  ENDIF()

  #Simple-Templated
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Simple-Templated
    SOURCES Simple-Templated.cpp
    DEPLIBS muelu muelu-xpetra
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

  #Simple-Belos
  IF(${PACKAGE_NAME}_ENABLE_Belos)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Simple-Belos
      SOURCES Simple-Belos.cpp
      NUM_MPI_PROCS 4
      COMM serial mpi
      DEPLIBS muelu muelu-xpetra muelu-adapters
      )
  ENDIF()

  #Simple-AztecOO
  IF(Trilinos_ENABLE_AztecOO) 
    # Use Trilinos_ENABLE_AztecOO instead of MueLu_ENABLE_AztecOO because this example requires AztecOO but MueLu doesn't.
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Simple-AztecOO
      SOURCES Simple-AztecOO.cpp
      NUM_MPI_PROCS 4
      COMM serial mpi
      DEPLIBS muelu muelu-xpetra muelu-adapters
      )
  ENDIF()

ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)

  #Simple-Tpetra
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Simple-Tpetra
    SOURCES Simple-Tpetra.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
    DEPLIBS muelu muelu-xpetra muelu-adapters
    )

  #Simple-Complex
  IF(HAVE_MUELU_COMPLEX)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Simple-Complex
      SOURCES Simple-Complex.cpp
      DEPLIBS muelu muelu-xpetra
      NUM_MPI_PROCS 1
      COMM serial mpi
      )
  ENDIF()

ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Ifpack AND ${PACKAGE_NAME}_ENABLE_Amesos)

  #Simple-Epetra
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Simple-Epetra
    SOURCES Simple-Epetra.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
    DEPLIBS muelu muelu-xpetra muelu-adapters
    )

ENDIF()

#TODO: NUM_MPI_PROCS->4
