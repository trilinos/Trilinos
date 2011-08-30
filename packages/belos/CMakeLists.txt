INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Belos)

#
# B) Set up package-specific options
#

# Kokkos includes an implementation of the Tall Skinny QR (TSQR)
# factorization.  TSQR is used to implement an OrthoManager
# (Belos::TsqrOrthoManager) and a MatOrthoManager (which only uses
# TSQR in the Euclidean-inner-product case).  TSQR is enabled by
# default if enabled in Kokkos.  If Teuchos_ENABLE_COMPLEX is enabled,
# TSQR is built with complex arithmetic support.
ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_TSQR
  HAVE_BELOS_TSQR
  "Enable TSQR (Tall Skinny QR factorization) support in Belos."
  "${Kokkos_ENABLE_TSQR}" )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR
  BELOS_TEUCHOS_TIME_MONITOR
  "Build with Teuchos TimeMonitors enabled."
  "${Trilinos_ENABLE_TEUCHOS_TIME_MONITOR}" )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Tpetra_Timers
  HAVE_BELOS_TPETRA_TIMERS
  "Add timers to Belos::MultiVecTraits for Tpetra multivectors."
  NO )

#
# "Experimental" code lives in Belos, but is not quite ready for use
# by a general audience.  You have to enable the CMake option below in
# order to install experimental headers, include experimental source
# files in the libraries, and build experimental tests.
#
# It's up to Belos developers to define what headers and source files
# are experimental, and what tests depend on experimental code.  They
# should disable building and running tests that have a required
# dependency on experimental code.  Some tests may have an optional
# dependency on experimental code; test implementers can test for this
# via "#ifdef HAVE_BELOS_EXPERIMENTAL".
#
ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Experimental
  HAVE_BELOS_EXPERIMENTAL
  "Enables experimental features of Belos, such as new solvers which are not quite ready for release to a general audience."
  NO )


#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

IF (Belos_ENABLE_Epetra)
  ADD_SUBDIRECTORY(epetra)
ENDIF()

IF (Belos_ENABLE_Tpetra)
  ADD_SUBDIRECTORY(tpetra)
ENDIF()


#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
