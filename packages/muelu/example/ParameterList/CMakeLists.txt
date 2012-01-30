INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Executables
#

IF((${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2) OR 
    (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Ifpack  AND ${PACKAGE_NAME}_ENABLE_Amesos))
  
  #MLCompatible
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MLCompatible
    SOURCES MLCompatible.cpp
    DEPLIBS muelu muelu-xpetra
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(MLCompatible_cp
    SOURCE_FILES ml_ParameterList.xml
    )
  
ENDIF()
