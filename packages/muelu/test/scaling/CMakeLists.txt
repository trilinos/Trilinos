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

INCLUDE(TribitsAddExecutable)

#TODO: how to remove that?
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)

# ScalingTest is compile with long long int when Teuchos_ENABLE_LONG_LONG_INT=ON
# The name of the executable is changed to make sure that long int support is correctly
# set when we do scalability tests.
IF(Teuchos_ENABLE_LONG_LONG_INT)
  IF(NOT ${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION) #FIXME
    #FIXME no Tpetra instantiation for long long int
    #FIXME we need a HAVE_MUELU_LONG_LONG_INT option
    SET(MueLu_ScalingTestName ScalingTest-LongLongInt)
  ELSE()
    SET(MueLu_ScalingTestName ScalingTest-Int)
  ENDIF()
ELSE()
  SET(MueLu_ScalingTestName ScalingTest-Int)
ENDIF()
#TODO: It is actually possible to build and test ScalingTest-Int also when LONG_LONG_INT=ON

IF (${PACKAGE_NAME}_ENABLE_Galeri)
IF ((${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2) OR
    (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Ifpack  AND ${PACKAGE_NAME}_ENABLE_Amesos))

  TRIBITS_ADD_EXECUTABLE(
    ${MueLu_ScalingTestName}
    SOURCES ScalingTest.cpp
    COMM serial mpi
    DEPLIBS muelu xpetra muelu-adapters
    )

  TRIBITS_ADD_EXECUTABLE( #TODO: to be add to nightly tests
    ScalingTestParamList  #TODO: Suffix Int/LongInt when done in the source code
    SOURCES ScalingTestParamList.cpp
    COMM serial mpi
    DEPLIBS muelu xpetra muelu-adapters
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(ScalingTestParamList_cp
    SOURCE_FILES scalingTest.xml
    )

ENDIF()

ENDIF() # Galeri

#
# Tests
#

INCLUDE(XpetraAddTest)

IF (${PACKAGE_NAME}_ENABLE_Galeri)

IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Ifpack  AND ${PACKAGE_NAME}_ENABLE_Amesos AND (NOT Teuchos_ENABLE_LONG_LONG_INT))
  XPETRA_ADD_TEST_EPETRA(${MueLu_ScalingTestName} 4)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)
#  XPETRA_ADD_TEST_TPETRA(${MueLu_ScalingTestName} 4)
ENDIF()

ENDIF() # Galeri