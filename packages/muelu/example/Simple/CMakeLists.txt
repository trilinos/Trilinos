INCLUDE(PackageAddExecutable)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Executables
#

IF(${PACKAGE_NAME}_ENABLE_Tpetra OR ${PACKAGE_NAME}_ENABLE_Epetra)

  #Simple
  PACKAGE_ADD_EXECUTABLE(
    Simple
    SOURCES Simple.cpp
    DEPLIBS muelu muelu-xpetra
    )

  #Simple-LongInt
  IF(Teuchos_ENABLE_LONG_LONG_INT)
    PACKAGE_ADD_EXECUTABLE(
      Simple-LongInt
      SOURCES Simple-LongInt.cpp
      DEPLIBS muelu muelu-xpetra
      )
  ENDIF()

  #Simple-Templated
  PACKAGE_ADD_EXECUTABLE(
    Simple-Templated
    SOURCES Simple-Templated.cpp
    DEPLIBS muelu muelu-xpetra
    )

  #Simple-Belos
  IF(${PACKAGE_NAME}_ENABLE_Belos)
    PACKAGE_ADD_EXECUTABLE(
      Simple-Belos
      SOURCES Simple-Belos.cpp
      DEPLIBS muelu muelu-xpetra muelu-adapters
      )
  ENDIF()

ENDIF()

#
# Tests
#

INCLUDE(../../xpetra/cmake/AddTestTpetraAndEpetra.cmake)

#Simple
PACKAGE_ADD_TEST_SEQ_TPETRA_AND_EPETRA(Simple)

#Simple-LongInt
IF(Teuchos_ENABLE_LONG_LONG_INT)
  PACKAGE_ADD_TEST_SEQ_TPETRA_AND_EPETRA(Simple-LongInt) 
ENDIF()

#Simple-Templated
PACKAGE_ADD_TEST_SEQ_TPETRA_AND_EPETRA(Simple-Templated)

#Simple-Belos
PACKAGE_ADD_TEST_SEQ_TPETRA_AND_EPETRA(Simple-Belos)

