INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(XpetraAddTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF ({PACKAGE_NAME}_ENABLE_Galeri)

#
#
#

IF((${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2) OR 
    (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Ifpack AND ${PACKAGE_NAME}_ENABLE_Amesos))

  #
  # ParameterList
  #

  TRIBITS_ADD_EXECUTABLE(
    ParameterList
    SOURCES ParameterList.cpp
    DEPLIBS muelu muelu-xpetra
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(ParameterList_cp
    SOURCE_FILES muelu_ParameterList.xml
    )

  #
  # MLParameterList  
  #

  TRIBITS_ADD_EXECUTABLE(
    MLParameterList
    SOURCES MLParameterList.cpp
    DEPLIBS muelu muelu-xpetra
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(MLParameterList_cp
    SOURCE_FILES ml_ParameterList.xml
    )
    
ENDIF()

#
# Tests of ParameterList and MLParameterList
#

IF(${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)
  XPETRA_ADD_TEST_TPETRA(ParameterList   1-4)
  XPETRA_ADD_TEST_TPETRA(MLParameterList 1-4)
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Ifpack AND ${PACKAGE_NAME}_ENABLE_Amesos)
  XPETRA_ADD_TEST_EPETRA(ParameterList   1-4)
  XPETRA_ADD_TEST_EPETRA(MLParameterList 1-4)
ENDIF()

#
# Stratimikos
#

IF((${PACKAGE_NAME}_ENABLE_Stratimikos) AND
    (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Ifpack AND ${PACKAGE_NAME}_ENABLE_Amesos))
  
  TRIBITS_ADD_EXECUTABLE(
    Stratimikos
    SOURCES Stratimikos.cpp
    DEPLIBS muelu muelu-xpetra muelu-adapters
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Stratimikos_cp
    SOURCE_FILES stratimikos_ParameterList.xml stratimikos_ParameterList2.xml
    )

  TRIBITS_ADD_TEST(
    Stratimikos
    NAME Stratimikos
    NUM_MPI_PROCS 4
    COMM serial mpi
    )
  
  TRIBITS_ADD_TEST(
    Stratimikos
    NAME Stratimikos2
    ARGS "--xml=stratimikos_ParameterList2.xml"
    NUM_MPI_PROCS 4
    COMM serial mpi
    )
  
ENDIF()

ENDIF()