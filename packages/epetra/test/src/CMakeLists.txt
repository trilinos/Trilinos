INCLUDE(TribitsLibraryMacros)

INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/test/src)

SET(TEST_HEADERS
  ${PACKAGE_SOURCE_DIR}/test/src/Epetra_test_functions.h
  ${PACKAGE_SOURCE_DIR}/test/src/Epetra_matrix_data.h
)

SET(TEST_SOURCES
  ${PACKAGE_SOURCE_DIR}/test/src/Epetra_test_functions.cpp
  ${PACKAGE_SOURCE_DIR}/test/src/Epetra_matrix_data.cpp
)

TRIBITS_ADD_LIBRARY(
  epetratest
  HEADERS ${TEST_HEADERS}
  SOURCES ${TEST_SOURCES}
  DEFINES -DEPETRA_LIB_EXPORTS_MODE
  DEPLIBS epetra
  TESTONLY
  )
