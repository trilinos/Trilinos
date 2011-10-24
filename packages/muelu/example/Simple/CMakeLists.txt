INCLUDE(PackageAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(${PACKAGE_NAME}_ENABLE_Tpetra OR ${PACKAGE_NAME}_ENABLE_Epetra)

  PACKAGE_ADD_EXECUTABLE(
    Simple
    SOURCES Simple.cpp
    COMM serial mpi
    DEPLIBS muelu muelu-xpetra
    )

  IF(Teuchos_ENABLE_LONG_LONG_INT)

    PACKAGE_ADD_EXECUTABLE(
      Simple-LongInt
      SOURCES Simple-LongInt.cpp
      COMM serial mpi
      DEPLIBS muelu muelu-xpetra
      )
    
  ENDIF()

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Tpetra)

  PACKAGE_ADD_TEST(
    Simple
    NAME "SimpleTpetra"
    ARGS "--linAlgebra=1"
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Epetra)

  PACKAGE_ADD_TEST(
    Simple
    NAME "SimpleEpetra"
    ARGS "--linAlgebra=0"
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

ENDIF()
