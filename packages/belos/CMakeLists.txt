INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

######################################################################
## Define the package
######################################################################

PACKAGE(Belos)

######################################################################
## Set up package-specific options
######################################################################

# Kokkos includes an implementation of the Tall Skinny QR (TSQR)
# factorization.  TSQR is used to implement orthogonalization methods.
# TSQR is enabled by default in Belos if enabled in Kokkos.  If
# Teuchos_ENABLE_COMPLEX is enabled, TSQR is built with complex
# arithmetic support.
ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_TSQR
  HAVE_BELOS_TSQR
  "Enable TSQR (Tall Skinny QR factorization) support in Belos."
  "${Kokkos_ENABLE_TSQR}" )

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_TEUCHOS_TIME_MONITOR
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
# Belos developers must maintain the list of experimental header and
# source files in this file, in Belos_EXPERIMENTAL_HEADERS
# resp. Belos_EXPERIMENTAL_SOURCES.  (This is a bit awkward because
# those files live in belos/src and should really be governed by
# belos/src/CMakeLists.txt.  However, this file needs to know about
# them so that it can exclude them from the release tarball.
# belos/src/CMakeLists.txt will reference those variables in order to
# exclude experimental code when building Belos' library/ies.)
#
# Belos developers must also exclude directories containing
# experimental tests (i.e., tests with required dependencies on
# experimental code) from both the release tarball (in the
# PACKAGE_ARCH_EXCLUDE_FILES section below) and the build (in the
# CMakeLists.txt files in the parent directory of the tests in
# question).  For the latter, look at tpetra/test/CMakeLists.txt for
# an example.
#
ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Experimental
  HAVE_BELOS_EXPERIMENTAL
  "Enables experimental features of Belos, such as new solvers which are not quite ready for release to a general audience."
  NO )


######################################################################
## Add the libraries, tests, and examples
######################################################################

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

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
PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

#
# List of all header files in belos/src that are marked
# "experimental."  See the comments above on "experimental" code.
#
SET(${PACKAGE_NAME}_EXPERIMENTAL_HEADERS "")
APPEND_SET(${PACKAGE_NAME}_EXPERIMENTAL_HEADERS
  BelosAkxFactory.hpp
  BelosAkx.hpp
  BelosBlockGCRODRIter.hpp
  BelosBlockGCRODRSolMgr.hpp
  BelosCaGmres.hpp
  BelosGmresBaseFactory.hpp
  BelosGmresBase.hpp
  BelosGmresBaseIteration.hpp
  BelosGmresInnerSolver.hpp
  BelosGmresSolMgr.hpp
  BelosInnerSolveResult.hpp
  BelosInnerSolver.hpp
  BelosMonomialOpAkx.hpp
  BelosNewtonBasis.hpp
  BelosNewtonOpAkx.hpp
  BelosOpAkx.hpp
  BelosStandardGmres.hpp
)

#
# List of all source files in belos/src that are marked
# "experimental."
#
SET(${PACKAGE_NAME}_EXPERIMENTAL_SOURCES "")
APPEND_SET(${PACKAGE_NAME}_EXPERIMENTAL_SOURCES
  BelosInnerSolveResult.cpp
)

#
# Exclude experimental header and source files, along with test
# directories for experimental code (i.e., tests with required
# dependencies on experimental code), from the release tarball.
#
PACKAGE_ARCH_EXCLUDE_FILES(
  ${${PACKAGE_NAME}_EXPERIMENTAL_HEADERS}
  ${${PACKAGE_NAME}_EXPERIMENTAL_SOURCES}
  tpetra/test/Gmres
  )

######################################################################
## Do standard postprocessing
######################################################################

PACKAGE_POSTPROCESS()
