INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../utils)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)

# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Epetra
  ${PACKAGE_NAME}_ENABLE_Tpetra
  )

# ScalingTest is compile with long long int when Teuchos_ENABLE_LONG_LONG_INT=ON
# The name of the executable is changed to make sure that long int support is correctly 
# set when we do scalability tests.
IF(Teuchos_ENABLE_LONG_LONG_INT)
  SET(MueLu_ScalingTestName ScalingTest-LongLongInt)
ELSE()
  SET(MueLu_ScalingTestName ScalingTest-Int)
ENDIF()
#TODO: It is actually possible to build ScalingTest-Int in both case

#
# Executable
#

PACKAGE_ADD_EXECUTABLE(
  ${MueLu_ScalingTestName}
  SOURCES ScalingTest.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-xpetra muelu-adapters
  )

#
# Tests
#

IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Amesos)

  PACKAGE_ADD_TEST(
    ${MueLu_ScalingTestName}
    NAME "ScalingTestEpetra"
    ARGS "--linAlgebra=0"
    NUM_MPI_PROCS 1
    COMM serial mpi
    )

  PACKAGE_ADD_TEST(
    ${MueLu_ScalingTestName}
    NAME "ScalingTestEpetra"
    ARGS "--linAlgebra=0"
    NUM_MPI_PROCS 4
    COMM mpi
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Amesos2)
  
  PACKAGE_ADD_TEST(
    ${MueLu_ScalingTestName}
    NAME "ScalingTestTpetra"
    ARGS "--linAlgebra=1"
    NUM_MPI_PROCS 1
    COMM serial mpi
    )

  PACKAGE_ADD_TEST(
    ${MueLu_ScalingTestName}
    NAME "ScalingTestTpetra"
    ARGS "--linAlgebra=1"
    NUM_MPI_PROCS 4
    COMM mpi
    )

ENDIF()

