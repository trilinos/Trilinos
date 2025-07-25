
INCLUDE(TrilinosCreateClientTemplateHeaders)
INCLUDE(ETI_functions)

#
# A) Package-specific configuration options
#

# TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

# SET(HEADERS ${HEADERS}
#   ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
#   )

IF (${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
  SET (MultiVector_ETI_SCALARS ${TpetraCore_ETI_SCALARS})
  IF (NOT Tpetra_INST_INT_INT)
    # GO = int is disabled, so the Scalar = GlobalOrdinal case doesn't
    # cover Scalar = int.  We need Scalar = int for communication of
    # process ranks.
    LIST(APPEND MultiVector_ETI_SCALARS "int")
  ENDIF()

  # Tpetra::CrsMatrix ETI Support
  ASSERT_DEFINED(Tpetra_ENABLE_Reduced_ETI)
  IF(Tpetra_ENABLE_Reduced_ETI)
    # For Reduced_ETI we don't include int
    SET (CrsMatrix_ETI_SCALARS  ${TpetraCore_ETI_SCALARS})
  ELSE()
    # Zoltan2 wants Scalar = int for CrsMatrix (Bug 6298), and thus
    # indirectly for RowMatrix as well.
    SET (CrsMatrix_ETI_SCALARS ${MultiVector_ETI_SCALARS})
  ENDIF()

  # Generate ETI .cpp files for Tpetra::MatrixMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(MATRIXMATRIX_OUTPUT_FILES "TpetraExt_ETI_SC_LO_GO_NT.tmpl" "MatrixMatrix" "MATRIXMATRIX" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${MATRIXMATRIX_OUTPUT_FILES})

ENDIF()

#
# Core Files
#

# general headers, sources
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# Must glob the binary dir last to get all of the auto-generated headers
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  tpetraext
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS tpetra
  )

# TRIBITS_ADD_TEST_DIRECTORIES(test)
