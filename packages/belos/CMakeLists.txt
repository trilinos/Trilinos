INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

######################################################################
## Define the package
######################################################################

TRIBITS_PACKAGE(Belos)

######################################################################
## Set up package-specific options
######################################################################

# Kokkos includes an implementation of the Tall Skinny QR (TSQR)
# factorization.  TSQR is used to implement orthogonalization methods.
# TSQR is enabled by default in Belos if enabled in Kokkos.  If
# Teuchos_ENABLE_COMPLEX is enabled, TSQR is built with complex
# arithmetic support.
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_TSQR
  HAVE_BELOS_TSQR
  "Enable TSQR (Tall Skinny QR factorization) support in Belos."
  "${KokkosClassic_ENABLE_TSQR}" )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR
  BELOS_TEUCHOS_TIME_MONITOR
  "Build with Teuchos TimeMonitors enabled."
  "${${PROJECT_NAME}_ENABLE_TEUCHOS_TIME_MONITOR}" )

TRIBITS_ADD_OPTION_AND_DEFINE(
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
# must use ${PACKAGE_NAME}_ENABLE_Experimental to conditionally
# disable tests or examples that have a required dependency on
# experimental code.  (${PACKAGE_NAME} evaluates to Belos, but we
# write it out so you know what to look for in the file.)  Code that
# _optionally_ depends on experimental code must be protected by
# #ifdef HAVE_BELOS_EXPERIMENTAL ... #endif.
#
# Belos developers must maintain the lists of experimental header
# resp. source files in belos/src/CMakeLists.txt, as
# ${PACKAGE_NAME}_EXPERIMENTAL_HEADERS
# resp. ${PACKAGE_NAME}_EXPERIMENTAL_SOURCES.
#
# Note that some experimental files may be excluded from the release
# tarball.  Belos uses the TRIBITS_EXCLUDE_FILES macro later in
# this file to list the files and directories to exclude from the
# release tarball.  That macro must NOT include conditional
# expressions; it must evaluate to the same list, no matter how
# Trilinos is configured.  We ensure this in Belos by listing
# explicitly the files and directories to exclude from the release
# tarball, rather than using variables or conditional CMake
# expressions.
#
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Experimental
  HAVE_BELOS_EXPERIMENTAL
  "Enable building and installation of experimental Belos features."
  NO )

# See above discussion of "experimental" Belos headers and source
# files.  belos/src/CMakeLists.txt may append to these lists.
SET(${PACKAGE_NAME}_EXPERIMENTAL_HEADERS "")
SET(${PACKAGE_NAME}_EXPERIMENTAL_SOURCES "")


######################################################################
## Add the libraries, tests, and examples
######################################################################

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

IF (Belos_ENABLE_Epetra)
  ADD_SUBDIRECTORY(epetra)
ENDIF()

IF (Belos_ENABLE_Tpetra)
  ADD_SUBDIRECTORY(tpetra)
ENDIF()


######################################################################
## Build list of files to exclude from the release tarball
######################################################################

#
# Exclude leftover GNU Autotools files.
#
TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# Files and directories to exclude from the release tarball.
#
TRIBITS_EXCLUDE_FILES(
  src/BelosAkxFactory.hpp
  src/BelosAkx.hpp
  src/BelosCaGmres.hpp
  src/BelosGmresBaseFactory.hpp
  src/BelosGmresBase.hpp
  src/BelosGmresBaseIteration.hpp
  src/BelosGmresInnerSolver.hpp
  src/BelosGmresSolMgr.hpp
  src/BelosInnerSolveResult.hpp
  src/BelosInnerSolver.hpp
  src/BelosMonomialOpAkx.hpp
  src/BelosNewtonBasis.hpp
  src/BelosNewtonOpAkx.hpp
  src/BelosOpAkx.hpp
  src/BelosProjectedLeastSquaresSolver.hpp
  src/BelosStandardGmres.hpp
  src/BelosInnerSolveResult.cpp
  tpetra/test/Gmres
  test/ProjectedLeastSquaresSolver
)

######################################################################
## Do standard postprocessing
######################################################################

TRIBITS_PACKAGE_POSTPROCESS()
