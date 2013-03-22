
INCLUDE(TribitsLibraryMacros)

######################################################################
## Package-specific configuration options
######################################################################

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

######################################################################
## Define the header and source files (and directories)
######################################################################

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# CORE
APPEND_SET(HEADERS
  BelosBlockCGIter.hpp
  BelosBlockCGSolMgr.hpp
  BelosBlockFGmresIter.hpp
  BelosBlockGmresIter.hpp
  BelosBlockGmresSolMgr.hpp
  BelosCGIteration.hpp
  BelosCGIter.hpp
  BelosConfigDefs.hpp
  BelosDGKSOrthoManager.hpp
  BelosGCRODRIter.hpp
  BelosGCRODRSolMgr.hpp
  BelosGmresIteration.hpp
  BelosGmresPolyOp.hpp
  BelosGmresPolySolMgr.hpp
  BelosICGSOrthoManager.hpp
  BelosIMGSOrthoManager.hpp
  BelosIteration.hpp
  BelosLinearProblem.hpp
  BelosLSQRIter.hpp
  BelosLSQRIteration.hpp
  BelosLSQRSolMgr.hpp
  BelosLSQRStatusTest.hpp
  BelosMatOrthoManager.hpp
  BelosMinresIteration.hpp
  BelosMinresIter.hpp
  BelosMinresSolMgr.hpp
  BelosMultiVec.hpp
  BelosMultiVecTraits.hpp
  BelosMVOPTester.hpp
  BelosOperator.hpp
  BelosOperatorTraits.hpp
  BelosOrthoManagerFactory.hpp
  BelosOrthoManager.hpp
  BelosOrthoManagerTest.hpp
  BelosOutputManager.hpp
  BelosPCPGIter.hpp
  BelosPCPGSolMgr.hpp
  BelosPseudoBlockCGIter.hpp
  BelosPseudoBlockCGSolMgr.hpp
  BelosPseudoBlockGmresIter.hpp
  BelosPseudoBlockGmresSolMgr.hpp
  BelosPseudoBlockStochasticCGIter.hpp
  BelosPseudoBlockStochasticCGSolMgr.hpp
  BelosRCGIter.hpp
  BelosRCGSolMgr.hpp
  BelosSimpleOrthoManager.hpp
  BelosSolverManager.hpp
  BelosSolverFactory.hpp
  BelosStatusTestCombo.hpp
  BelosStatusTestFactory.hpp
  BelosStatusTestGeneralOutput.hpp
  BelosStatusTestGenResNorm.hpp
  BelosStatusTest.hpp
  BelosStatusTestImpResNorm.hpp
  BelosStatusTestMaxIters.hpp
  BelosStatusTestOutputFactory.hpp
  BelosStatusTestOutput.hpp
  BelosStatusTestResNorm.hpp
  BelosStatusTestResNormOutput.hpp
  BelosStochasticCGIteration.hpp
  BelosStubTsqrAdapter.hpp
  BelosTFQMRIter.hpp
  BelosTFQMRSolMgr.hpp
  BelosTypes.hpp
  ) 

# UTIL
APPEND_SET(SOURCES
  BelosTypes.cpp
  BelosVersion.cpp
  )

# COMP_H
APPEND_SET(HEADERS
  BelosBlockCGSolMgr.hpp
  BelosMinresSolMgr.hpp
  BelosBlockGmresSolMgr.hpp
  BelosPseudoBlockCGSolMgr.hpp
  BelosPseudoBlockGmresSolMgr.hpp
  BelosPseudoBlockStochasticCGSolMgr.hpp
  BelosTFQMRSolMgr.hpp
  BelosCGIter.hpp
  BelosBlockCGIter.hpp
  BelosMinresIter.hpp
  BelosBlockGmresIter.hpp
  BelosPseudoBlockCGIter.hpp
  BelosPseudoBlockGmresIter.hpp
  BelosPseudoBlockStochasticCGIter.hpp
  BelosTFQMRIter.hpp
  )

######################################################################
## Add files to HEADERS and /or SOURCES, depending on the values of
## certain configure-time options.
######################################################################

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_TSQR)
IF(${PACKAGE_NAME}_ENABLE_TSQR)
  #
  # Build the TSQR-based OrthoManager and MatOrthoManager classes.
  #
  APPEND_SET(HEADERS BelosTsqrOrthoManagerImpl.hpp BelosTsqrOrthoManager.hpp)
ENDIF() # ${PACKAGE_NAME}_ENABLE_TSQR

#
# See the comments on "experimental" code in belos/CMakeLists.txt.  In
# particular, experimental headers and source files are not included
# in the build or install unless the Belos_ENABLE_Experimental CMake
# option was set at configure time.
#
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
  BelosProjectedLeastSquaresSolver.hpp
  BelosStandardGmres.hpp
)

APPEND_SET(${PACKAGE_NAME}_EXPERIMENTAL_SOURCES
  BelosInnerSolveResult.cpp
)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Experimental)
IF(${PACKAGE_NAME}_ENABLE_Experimental)
  #
  # Add experimental files to HEADERS and SOURCES.
  # 
  APPEND_SET(HEADERS ${${PACKAGE_NAME}_EXPERIMENTAL_HEADERS})
  APPEND_SET(SOURCES ${${PACKAGE_NAME}_EXPERIMENTAL_SOURCES})
ENDIF() # ${PACKAGE_NAME}_ENABLE_Experimental

######################################################################
## Define the targets for the package's library/ies.
######################################################################

IF (${PACKAGE_NAME}_ENABLE_TSQR)
  # If TSQR support is enabled in Belos, we have to link Belos with
  # the appropriate Kokkos libraries, since TSQR is implemented in
  # Kokkos.
  TRIBITS_ADD_LIBRARY(
    belos
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS kokkos kokkosnodetsqr kokkosdisttsqr
    )
ELSE ()
  TRIBITS_ADD_LIBRARY(
    belos
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )
ENDIF ()
