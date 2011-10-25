INCLUDE(PackageAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Executables
#

#TODO: - NUM_MPI_PROCS == 1 because examples are using a sequential smoother
#      - Test for Belos disabled

IF(${PACKAGE_NAME}_ENABLE_Tpetra OR ${PACKAGE_NAME}_ENABLE_Epetra)

  #Simple
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Simple
    SOURCES Simple.cpp
    DEPLIBS muelu muelu-xpetra
    NUM_MPI_PROCS 1
    COMM serial mpi
    )

  #Simple-LongInt
  IF(Teuchos_ENABLE_LONG_LONG_INT)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      Simple-LongInt
      SOURCES Simple-LongInt.cpp
      DEPLIBS muelu muelu-xpetra
      NUM_MPI_PROCS 1
      COMM serial mpi
      )
  ENDIF()

  #Simple-Templated
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Simple-Templated
    SOURCES Simple-Templated.cpp
    DEPLIBS muelu muelu-xpetra
    NUM_MPI_PROCS 1
    COMM serial mpi
    )

  IF(${PACKAGE_NAME}_ENABLE_Tpetra) #TMP
  #Simple-Belos
  IF(${PACKAGE_NAME}_ENABLE_Belos)
    PACKAGE_ADD_EXECUTABLE(
      Simple-Belos
      SOURCES Simple-Belos.cpp
      DEPLIBS muelu muelu-xpetra muelu-adapters
      )
  ENDIF()
  ENDIF()

ENDIF()
