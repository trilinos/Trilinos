INCLUDE(TribitsAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Executables
#

IF((${PACKAGE_NAME}_ENABLE_Tpetra  AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2) OR
    (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND ${PACKAGE_NAME}_ENABLE_Ifpack  AND ${PACKAGE_NAME}_ENABLE_Amesos))

  IF(${PACKAGE_NAME}_ENABLE_Galeri)

    TRIBITS_ADD_EXECUTABLE(
      FixedMatrixPattern
      SOURCES FixedMatrixPattern.cpp
      DEPLIBS muelu xpetra
      COMM serial mpi
      )

    TRIBITS_ADD_EXECUTABLE(
      StandardReuse
      SOURCES StandardReuse.cpp
      DEPLIBS muelu xpetra
      COMM serial mpi
      )
  ENDIF()

ENDIF()

#
# Tests
#

INCLUDE(XpetraAddTest)

IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2 AND ${PACKAGE_NAME}_ENABLE_Galeri)
  XPETRA_ADD_TEST_TPETRA(FixedMatrixPattern 4)
  XPETRA_ADD_TEST_TPETRA(StandardReuse 4)
ENDIF()
IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Ifpack  AND ${PACKAGE_NAME}_ENABLE_Amesos AND ${PACKAGE_NAME}_ENABLE_Galeri)
  XPETRA_ADD_TEST_EPETRA(FixedMatrixPattern 4)
  XPETRA_ADD_TEST_EPETRA(StandardReuse 4)
ENDIF()
