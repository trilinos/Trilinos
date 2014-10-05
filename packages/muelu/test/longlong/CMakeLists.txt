# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Epetra
  ${PACKAGE_NAME}_ENABLE_Tpetra
  ${PACKAGE_NAME}_ENABLE_Zoltan2
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
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

INCLUDE(XpetraAddTest)

IF(${PACKAGE_NAME}_ENABLE_Galeri AND ${PACKAGE_NAME}_ENABLE_Tpetra AND
   ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2 AND
   ${PACKAGE_NAME}_ENABLE_Zoltan2 AND Zoltan2_ENABLE_Experimental AND
   TPL_ENABLE_MPI)
  TRIBITS_ADD_EXECUTABLE(
    TwoBillionTest
    SOURCES TwoBillion.cpp
    COMM mpi
    DEPLIBS muelu xpetra muelu-adapters
    )


  TRIBITS_ADD_TEST(
    TwoBillionTest
    NAME "TwoBillionTest"
    NUM_MPI_PROCS 4
    COMM mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(TwoBillion_cp
    MESSAGE("doing copy")
    SOURCE_FILES TwoBillion.xml
    )


ENDIF()

