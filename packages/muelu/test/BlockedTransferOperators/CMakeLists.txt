INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)

IF ((${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 ) OR 
    (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Ifpack  ))

	TRIBITS_ADD_EXECUTABLE_AND_TEST(
	    BlockedTransferOperators
	    SOURCES BlockedTransferOperators.cpp
	    ARGS ""
	    COMM serial mpi
	    DEPLIBS muelu muelu-xpetra
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
#    DEPLIBS muelu muelu-xpetra
#  )
  
#TRIBITS_COPY_FILES_TO_BINARY_DIR(TSIExample_test_copyfiles
#SOURCE_FILES TSImatrix.out TSIbvector.out nspvector.out
#)
  
#ENDIF()