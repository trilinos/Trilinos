
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE( Anasazi ENABLE_SHADOWING_WARNINGS )

# 2009/05/01: rabartl: I have removed 'CLEANED' from above from now because it
# fails to compile with gcc 3.4.6 (see bug 4395).

#
# B) Set up package-specific options
#

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

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_COMPLEX
  HAVE_ANASAZI_COMPLEX
  "Enable Anasazi features for complex data types."
  "${Teuchos_ENABLE_COMPLEX}" )

# Kokkos includes an implementation of the Tall Skinny QR (TSQR)
# factorization.  TSQR is used to implement an OrthoManager
# (Anasazi::TsqrOrthoManager) and a MatOrthoManager (which only uses
# TSQR in the Euclidean-inner-product case).  TSQR is enabled by
# default if enabled in Kokkos.  If Teuchos_ENABLE_COMPLEX is enabled,
# TSQR is built with complex arithmetic support.

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_TSQR
  HAVE_ANASAZI_TSQR
  "Enable TSQR (Tall Skinny QR factorization) support in Anasazi."
  "${Kokkos_ENABLE_TSQR}" )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR
  ANASAZI_TEUCHOS_TIME_MONITOR
  "Build with Teuchos TimeMonitors enabled."
  "${Trilinos_ENABLE_TEUCHOS_TIME_MONITOR}" )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Tpetra_Timers
  HAVE_ANASAZI_TPETRA_TIMERS
  "Add timers to Anasazi::MultiVecTraits for Tpetra multivectors."
  NO
  )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

IF (Anasazi_ENABLE_Epetra)
  ADD_SUBDIRECTORY(epetra)
ENDIF()

PACKAGE_ADD_TEST_DIRECTORIES(test)

IF (Anasazi_ENABLE_ThyraCore)
  ADD_SUBDIRECTORY(thyra)
ENDIF()

IF (Anasazi_ENABLE_Tpetra)
  ADD_SUBDIRECTORY(tpetra)
ENDIF()

#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

PACKAGE_ARCH_EXCLUDE_FILES(
  doc/TOMS
  doc/OrthoStudy
  doc/ThyraPerf
  src/ModalAnalysisSolvers
  )

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
