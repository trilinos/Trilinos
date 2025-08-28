

#
# A) Define the package
#

TRIBITS_PACKAGE( Anasazi ENABLE_SHADOWING_WARNINGS )

# 2009/05/01: rabartl: I have removed 'CLEANED' from above from now because it
# fails to compile with gcc 3.4.6 (see bug 4395).

#
# B) Set up package-specific options
#

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

# RBGen support
TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_RBGen
  HAVE_ANASAZI_RBGEN
  "Enable Reduced Basis Generation (RBGen) support in Anasazi"
  OFF )

# Anasazi includes tests that may require complex arithmetic in order
# to build.  These include especially (but not exclusively) the
# Tpetra-specific tests, since Tpetra objects are templated on the
# Scalar data type.  We don't want to build tests that depend on
# complex arithmetic if the rest of Trilinos (esp. Teuchos, since we
# need e.g., Teuchos::ScalarTraits< Scalar >) was not built with
# complex arithmetic support.  Thus, we define an appropriate package
# configuration variable ${PACKAGE_NAME}_ENABLE_COMPLEX, which we can
# use to decide whether or not to build particular tests.
#
# The default value of this config variable should be the same as
# Teuchos_ENABLE_COMPLEX, because most of Anasazi's complex-arithmetic
# tests depend heavily on Teuchos' complex arithmetic support.

IF (Trilinos_ENABLE_COMPLEX_DOUBLE  AND  Teuchos_ENABLE_COMPLEX)
  SET(${PACKAGE_NAME}_ENABLE_COMPLEX_DEFAULT  ON)
ELSE()
  SET(${PACKAGE_NAME}_ENABLE_COMPLEX_DEFAULT OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_COMPLEX
  HAVE_ANASAZI_COMPLEX
  "Enable Anasazi features for complex data types."
  ${${PACKAGE_NAME}_ENABLE_COMPLEX_DEFAULT} )

# ToDo: Above, split ${PACKAGE_NAME}_ENABLE_COMPLEX into
# ${PACKAGE_NAME}_ENABLE_COMPLEX_FLOAT and
# ${PACKAGE_NAME}_ENABLE_COMPLEX_DOUBLE.  See Trilinos GitHub #362.

IF (${PACKAGE_NAME}_ENABLE_COMPLEX AND NOT Teuchos_ENABLE_COMPLEX)
  MESSAGE(FATAL_ERROR
    "Error: ${PACKAGE_NAME}_ENABLE_COMPLEX='${${PACKAGE_NAME}_ENABLE_COMPLEX}' "
    " and Teuchos_ENABLE_COMPLEX='${Teuchos_ENABLE_COMPLEX}'.  This would result in"
    " a build or link-time failure so we have to abort.")
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_COMPLEX AND
  ${PACKAGE_NAME}_ENABLE_Tpetra AND NOT Tpetra_INST_COMPLEX_DOUBLE
  )
  MESSAGE(FATAL_ERROR
    "Error: ${PACKAGE_NAME}_ENABLE_COMPLEX='${${PACKAGE_NAME}_ENABLE_COMPLEX}',"
    " ${PACKAGE_NAME}_ENABLE_Tpetra='${${PACKAGE_NAME}_ENABLE_Tpetra}',"
    " and Tpetra_INST_COMPLEX_DOUBLE='${Tpetra_INST_COMPLEX_DOUBLE}' ."
    "  This would result in a build or link-time failure so we have to abort.")
ENDIF()

# Kokkos includes an implementation of the Tall Skinny QR (TSQR)
# factorization.  TSQR is used to implement an OrthoManager
# (Anasazi::TsqrOrthoManager) and a MatOrthoManager (which only uses
# TSQR in the Euclidean-inner-product case).  TSQR is enabled by
# default if the KokkosTSQR and Tpetra packages are enabled, but users
# have control over Tpetra_ENABLE_TSQR.  Note that the Epetra adapter
# for TSQR currently lives in Tpetra.

IF (${PACKAGE_NAME}_ENABLE_Tpetra AND Tpetra_ENABLE_TSQR)
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_TSQR_DEFAULT ON)
ELSE()
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_TSQR_DEFAULT OFF)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_TSQR
  HAVE_ANASAZI_TSQR
  "Enable TSQR (Tall Skinny QR factorization) support in Anasazi."
  "${${PACKAGE_NAME}_ENABLE_TSQR_DEFAULT}")

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR
  ANASAZI_TEUCHOS_TIME_MONITOR
  "Build with Teuchos TimeMonitors enabled."
  "${${PROJECT_NAME}_ENABLE_TEUCHOS_TIME_MONITOR}" )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Tpetra_Timers
  HAVE_ANASAZI_TPETRA_TIMERS
  "Add timers to Anasazi::MultiVecTraits for Tpetra multivectors."
  NO
  )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

IF (${PACKAGE_NAME}_ENABLE_RBGen)
  ADD_SUBDIRECTORY(src-rbgen)
ENDIF()

IF (Anasazi_ENABLE_Epetra)
  ADD_SUBDIRECTORY(epetra)
ENDIF()

TRIBITS_ADD_TEST_DIRECTORIES(test)

IF (Anasazi_ENABLE_ThyraCore)
  ADD_SUBDIRECTORY(thyra)
ENDIF()

IF (Anasazi_ENABLE_Tpetra)
  ADD_SUBDIRECTORY(tpetra)
ENDIF()

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_FILES(
  epetra/src/ModalAnalysisSolvers
  )

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
