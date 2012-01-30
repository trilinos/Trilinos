INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(XpetraAddTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Executables
#

IF((${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2) OR 
    (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Ifpack AND ${PACKAGE_NAME}_ENABLE_Amesos))
  
  # MLCompatible
  TRIBITS_ADD_EXECUTABLE(
    MLCompatible
    SOURCES MLCompatible.cpp
    DEPLIBS muelu muelu-xpetra
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(MLCompatible_cp
    SOURCE_FILES ml_ParameterList.xml
    )
  
  IF (${PACKAGE_NAME}_ENABLE_Amesos2)
    XPETRA_ADD_TEST_TPETRA(MLCompatible "serial mpi" 1-4)
  ENDIF()
  
  IF (${PACKAGE_NAME}_ENABLE_Amesos)
    XPETRA_ADD_TEST_EPETRA(MLCompatible "serial mpi" 1-4)
  ENDIF()
  
ENDIF()
