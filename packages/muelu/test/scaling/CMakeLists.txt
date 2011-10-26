# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Epetra
  ${PACKAGE_NAME}_ENABLE_Tpetra
  )

#
# Executable
#

INCLUDE(PackageAddExecutable)

#TODO: how to remove that?
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../utils)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)

# ScalingTest is compile with long long int when Teuchos_ENABLE_LONG_LONG_INT=ON
# The name of the executable is changed to make sure that long int support is correctly 
# set when we do scalability tests.
IF(Teuchos_ENABLE_LONG_LONG_INT)
  SET(MueLu_ScalingTestName ScalingTest-LongLongInt)
ELSE()
  SET(MueLu_ScalingTestName ScalingTest-Int)
ENDIF()
#TODO: It is actually possible to build and test ScalingTest-Int also when LONG_LONG_INT=ON

PACKAGE_ADD_EXECUTABLE(
  ${MueLu_ScalingTestName}
  SOURCES ScalingTest.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-xpetra muelu-adapters
  )

#
# Tests
#

INCLUDE(XpetraAddTest)

IF (${PACKAGE_NAME}_ENABLE_Amesos)

  XPETRA_ADD_TEST_EPETRA(${MueLu_ScalingTestName} "mpi"        4)
  XPETRA_ADD_TEST_EPETRA(${MueLu_ScalingTestName} "serial mpi" 1)

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Amesos2)
  
  XPETRA_ADD_TEST_EPETRA(${MueLu_ScalingTestName} "mpi"        4)
  XPETRA_ADD_TEST_EPETRA(${MueLu_ScalingTestName} "serial mpi" 1)

ENDIF()
