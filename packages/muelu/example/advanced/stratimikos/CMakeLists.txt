INCLUDE(XpetraAddTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../test/unit_tests)

IF (${PACKAGE_NAME}_ENABLE_Galeri)

#
# Stratimikos example
#

IF(${PACKAGE_NAME}_ENABLE_Stratimikos)

# build executable
TRIBITS_ADD_EXECUTABLE(
  Stratimikos
  SOURCES Stratimikos.cpp
  COMM serial mpi
  )

IF(${PACKAGE_NAME}_ENABLE_Epetra AND
   ${PACKAGE_NAME}_ENABLE_Ifpack AND
   ${PACKAGE_NAME}_ENABLE_Amesos)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(advanced_stratimikos_epetra_cp
    SOURCE_FILES stratimikos_ParameterList.xml stratimikos_ParameterList2.xml
    )

  TRIBITS_ADD_TEST(
    Stratimikos
    NAME Stratimikos
    ARGS "--xml=stratimikos_ParameterList2.xml --linAlgebra=Epetra"
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
    Stratimikos
    NAME Stratimikos2
    ARGS "--xml=stratimikos_ParameterList2.xml --linAlgebra=Epetra"
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

ENDIF() # EPETRA solver stack

IF(${PACKAGE_NAME}_ENABLE_Tpetra AND
   ${PACKAGE_NAME}_ENABLE_Ifpack2 AND
   ${PACKAGE_NAME}_ENABLE_Amesos2)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(advanced_stratimikos_tpetra_cp
    SOURCE_FILES stratimikos_ParameterList.xml stratimikos_ParameterList2.xml
    )

  TRIBITS_ADD_TEST(
    Stratimikos
    NAME Stratimikos
    ARGS "--xml=stratimikos_ParameterList2.xml --linAlgebra=Tpetra"
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
    Stratimikos
    NAME Stratimikos2
    ARGS "--xml=stratimikos_ParameterList2.xml --linAlgebra=Tpetra"
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

ENDIF() # TPETRA solver stack

ENDIF() # STRATIMIKOS

ENDIF() # GALERI
