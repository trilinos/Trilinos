INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../utils)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)

# first ensure that these cmake boolean variables are defined
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos2)

PACKAGE_ADD_EXECUTABLE(
  ScalingTest
  SOURCES ScalingTest.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-xpetra muelu-adapters
  )

IF (${PACKAGE_NAME}_ENABLE_Epetra)
  IF (${PACKAGE_NAME}_ENABLE_Amesos)
    IF(NOT Teuchos_ENABLE_LONG_LONG_INT)

      PACKAGE_ADD_TEST(
        ScalingTest
        NAME "ScalingTestEpetra"
        ARGS "--linAlgebra=0"
        NUM_MPI_PROCS 1
        COMM serial mpi
        )

      PACKAGE_ADD_TEST(
        ScalingTest
        NAME "ScalingTestEpetra"
        ARGS "--linAlgebra=0"
        NUM_MPI_PROCS 4
        COMM mpi
        )
    ENDIF()

  ENDIF()
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Tpetra)
  IF (${PACKAGE_NAME}_ENABLE_Amesos2)
    
    PACKAGE_ADD_TEST(
      ScalingTest
      NAME "ScalingTestTpetra"
      ARGS "--linAlgebra=1"
      NUM_MPI_PROCS 1
      COMM serial mpi
      )

    PACKAGE_ADD_TEST(
      ScalingTest
      NAME "ScalingTestTpetra"
      ARGS "--linAlgebra=1"
      NUM_MPI_PROCS 4
      COMM mpi
      )

  ENDIF()
ENDIF()
