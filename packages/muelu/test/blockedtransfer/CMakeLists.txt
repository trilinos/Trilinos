TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)

TRIBITS_ADD_EXECUTABLE(
  BlockedTransfer
  SOURCES BlockedTransfer.cpp
  COMM serial mpi
)

IF(${PACKAGE_NAME}_ENABLE_Epetra
   AND ${PACKAGE_NAME}_ENABLE_EpetraExt
   AND ${PACKAGE_NAME}_ENABLE_Ifpack
)

  TRIBITS_ADD_TEST(
    BlockedTransfer
    NAME "BlockedTransfer_Epetra"
    ARGS "--linAlgebra=Epetra"
    COMM serial mpi
  )

ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Ifpack2)

  TRIBITS_ADD_TEST(
    BlockedTransfer
    NAME "BlockedTransfer_Tpetra"
    ARGS "--linAlgebra=Tpetra"
    COMM serial mpi
  )

ENDIF()

# first ensure that these cmake boolean variables are defined
#ASSERT_DEFINED(
#  ${PACKAGE_NAME}_ENABLE_Epetra
#  ${PACKAGE_NAME}_ENABLE_EpetraExt
#  ${PACKAGE_NAME}_ENABLE_Amesos
#  )

#IF (${PACKAGE_NAME}_ENABLE_Epetra AND
#    ${PACKAGE_NAME}_ENABLE_EpetraExt AND
#    ${PACKAGE_NAME}_ENABLE_Amesos)

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    TSIExample
#    SOURCES tsiexample.cpp
#    ARGS ""
#    COMM serial mpi
#  )

#TRIBITS_COPY_FILES_TO_BINARY_DIR(TSIExample_test_cp
#SOURCE_FILES TSImatrix.out TSIbvector.out nspvector.out
#)

#ENDIF()
