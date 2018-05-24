
######################################################################
## Define the package
######################################################################

TRIBITS_PACKAGE(Belos)

######################################################################
## Set up package-specific options
######################################################################

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

# Kokkos includes an implementation of the Tall Skinny QR (TSQR)
# factorization.  TSQR is used to implement an OrthoManager and a
# MatOrthoManager (which only uses TSQR in the Euclidean-inner-product
# case).  TSQR is enabled by default if the KokkosTSQR and Tpetra
# packages are enabled, but users have control over
# Tpetra_ENABLE_TSQR.  Note that the Epetra adapter for TSQR currently
# lives in Tpetra.
#
# mfh 05 Sep 2012: Back in 2010, I had to put Epetra's TSQR adapter in
# Tpetra.  This was because an external software package was linking
# to Epetra without using Trilinos' Makefile.export mechanism to
# consider which Trilinos libraries to use.  Making Epetra optionally
# depend on Kokkos (to make TSQR work) changed the list of libraries,
# and that software package wasn't doing the right thing to adjust.
# In any case, I chose at the time to put Epetra's TSQR adapter in
# Tpetra, since Tpetra already had an optional dependency on Epetra.
# So, for better or worse, you won't be able to use TSQR with Epetra
# unless you enable both Epetra and Tpetra in your Trilinos build.
# HAVE_BELOS_TSQR will correctly reflect this: it won't be defined 
# unless both Epetra and Tpetra are enabled in your Trilinos build.
IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Tpetra AND Tpetra_ENABLE_TSQR)
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_TSQR_DEFAULT ON)
ELSE()
  GLOBAL_SET(${PACKAGE_NAME}_ENABLE_TSQR_DEFAULT OFF)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_TSQR
  HAVE_BELOS_TSQR
  "Enable TSQR (Tall Skinny QR factorization) as an orthogonalization method in Belos.  You must build Trilinos with both Epetra _and_ Tpetra enabled in order to use TSQR with Epetra objects."
  "${${PACKAGE_NAME}_ENABLE_TSQR_DEFAULT}")

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR
  BELOS_TEUCHOS_TIME_MONITOR
  "Enable Teuchos timers (via TimeMonitor) in Belos solvers."
  ${${PROJECT_NAME}_ENABLE_TEUCHOS_TIME_MONITOR} )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Tpetra_Timers
  HAVE_BELOS_TPETRA_TIMERS
  "Add timers to Belos::MultiVecTraits for Tpetra::MultiVector operations."
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
  "Enable experimental Belos solvers and features."
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

IF (Belos_ENABLE_Epetra)
  ADD_SUBDIRECTORY(epetra)
ENDIF()

IF (Belos_ENABLE_Tpetra)
  ADD_SUBDIRECTORY(tpetra)
ENDIF()

IF (Belos_ENABLE_Xpetra)
  ADD_SUBDIRECTORY(xpetra)
ENDIF()

# Do this after building the libs - we need them
ADD_SUBDIRECTORY(doc)

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
  src/BelosInnerSolveResult.hpp
  src/BelosInnerSolver.hpp
  src/BelosProjectedLeastSquaresSolver.hpp
  src/BelosInnerSolveResult.cpp
  test/ProjectedLeastSquaresSolver
)

######################################################################
## Do standard postprocessing
######################################################################

TRIBITS_PACKAGE_POSTPROCESS()
