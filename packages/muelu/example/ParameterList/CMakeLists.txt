INCLUDE(XpetraAddTest)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
#
#

IF((${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)
   OR (${PACKAGE_NAME}_ENABLE_Epetra
       AND ${PACKAGE_NAME}_ENABLE_EpetraExt
       AND ${PACKAGE_NAME}_ENABLE_Ifpack
       AND ${PACKAGE_NAME}_ENABLE_Amesos
      )
)

  #
  # MLParameterList
  #

  TRIBITS_ADD_EXECUTABLE(
    MLParameterList
    SOURCES MLParameterList.cpp
    COMM serial mpi
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(MLParameterList_cp SOURCE_FILES ml_ParameterList.xml ml_RebalancedParameterList.xml)

ENDIF()

#
# Tests of MLParameterList
#

IF(${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)
  #  XPETRA_ADD_TEST_TPETRA(MLParameterList 1-4)
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Epetra
   AND ${PACKAGE_NAME}_ENABLE_EpetraExt
   AND ${PACKAGE_NAME}_ENABLE_Ifpack
   AND ${PACKAGE_NAME}_ENABLE_Amesos
)
  XPETRA_ADD_TEST_EPETRA(MLParameterList 1-4)

  TRIBITS_ADD_TEST(
    MLParameterList
    NAME MLParameterList_Repartition-Epetra
    ARGS "--xml=ml_RebalancedParameterList.xml --linAlgebra=Epetra --muelu2=0"
    NUM_MPI_PROCS 4
    COMM serial mpi
  )
ENDIF()

# TODO: add a test with the --xml option
