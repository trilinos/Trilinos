INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)

# first ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  )

ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos2
  )

IF (${PACKAGE_NAME}_ENABLE_Amesos2)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TransferOps
    SOURCES TransferOps.cpp
    ARGS ""
    COMM serial mpi
    DEPLIBS muelu xpetra muelu-adapters
  )

IF (${PACKAGE_NAME}_ENABLE_Epetra)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Condif2D
    SOURCES Condif2d.cpp
    ARGS ""
    COMM serial mpi
    DEPLIBS muelu xpetra muelu-adapters
  )


  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Condif2DFiltered
    SOURCES Condif2dFiltered.cpp
    ARGS ""
    COMM serial mpi
    DEPLIBS muelu xpetra muelu-adapters
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Condif2_cp
  SOURCE_FILES Condif2Mat.mat Condif2Rhs.mat
  )
ENDIF()

ENDIF()