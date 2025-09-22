TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../test/unit_tests)

#
# Executables
#

IF((${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)
   OR (${PACKAGE_NAME}_ENABLE_Epetra
       AND ${PACKAGE_NAME}_ENABLE_EpetraExt
       AND ${PACKAGE_NAME}_ENABLE_Ifpack
       AND ${PACKAGE_NAME}_ENABLE_Amesos
      )
)

  TRIBITS_ADD_EXECUTABLE(
    FixedMatrixPattern
    SOURCES FixedMatrixPattern.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    StandardReuse
    SOURCES StandardReuse.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    ReuseSequence
    SOURCES ReuseSequence.cpp
    COMM serial mpi
  )

ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(ReuseSequence_cp SOURCE_FILES reuse_seq.xml)

#
# Tests
#

INCLUDE(XpetraAddTest)

IF(${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)
  XPETRA_ADD_TEST_TPETRA(FixedMatrixPattern 4)
  XPETRA_ADD_TEST_TPETRA(StandardReuse 4)
  MUELU_ADD_SERIAL_AND_MPI_TEST(
    ReuseSequence
    NAME "ReuseSequenceTpetra"
    ARGS "--linAlgebra=Tpetra --xml=reuse_seq.xml --nx=40 --ny=40 --dim=2"
    NUM_MPI_PROCS 1
    COMM serial mpi
  )
ENDIF()
IF(${PACKAGE_NAME}_ENABLE_Epetra
   AND ${PACKAGE_NAME}_ENABLE_Ifpack
   AND ${PACKAGE_NAME}_ENABLE_Amesos
)
  XPETRA_ADD_TEST_EPETRA(FixedMatrixPattern 4)
  XPETRA_ADD_TEST_EPETRA(StandardReuse 4)
  MUELU_ADD_SERIAL_AND_MPI_TEST(
    ReuseSequence
    NAME "ReuseSequenceEpetra"
    ARGS "--linAlgebra=Epetra --xml=reuse_seq.xml --nx=40 --ny=40 --dim=2"
    NUM_MPI_PROCS 1
    COMM serial mpi
  )
ENDIF()
