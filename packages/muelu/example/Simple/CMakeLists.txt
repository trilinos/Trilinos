INCLUDE(PackageAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Executables
#

IF(${PACKAGE_NAME}_ENABLE_Tpetra OR ${PACKAGE_NAME}_ENABLE_Epetra)

  PACKAGE_ADD_EXECUTABLE(
    Simple
    SOURCES Simple.cpp
    DEPLIBS muelu muelu-xpetra
    )

  IF(Teuchos_ENABLE_LONG_LONG_INT)

    PACKAGE_ADD_EXECUTABLE(
      Simple-LongInt
      SOURCES Simple-LongInt.cpp
      DEPLIBS muelu muelu-xpetra
      )
    
  ENDIF()

ENDIF()

#
# Tests
#

INCLUDE(../../xpetra/cmake/AddTestTpetraAndEpetra.cmake)

PACKAGE_ADD_TEST_TPETRA_AND_EPETRA(Simple)
PACKAGE_ADD_TEST_TPETRA_AND_EPETRA(Simple-LongInt)
