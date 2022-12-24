

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
  Belos_Details_EBelosSolverType.hpp
  BelosBiCGStabIter.hpp
  BelosBiCGStabSolMgr.hpp
  BelosBlockCGIter.hpp
  BelosBlockCGSolMgr.hpp
  BelosBlockFGmresIter.hpp
  BelosBlockGmresIter.hpp
  BelosBlockGmresSolMgr.hpp
  BelosCGIteration.hpp
  BelosCGIter.hpp
  BelosCGSingleRedIter.hpp
  BelosConfigDefs.hpp
  BelosDGKSOrthoManager.hpp
  BelosFixedPointIter.hpp
  BelosFixedPointIteration.hpp
  BelosFixedPointSolMgr.hpp
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
  BelosOperatorT.hpp
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
  BelosSolverFactory_Belos.hpp
  BelosSolverFactory_Generic.hpp
  BelosCustomSolverFactory.hpp
  BelosStatusTestCombo.hpp
  BelosStatusTestFactory.hpp
  BelosStatusTestGeneralOutput.hpp
  BelosStatusTestGenResNorm.hpp
  BelosStatusTestGenResSubNorm.hpp
  BelosStatusTest.hpp
  BelosStatusTestImpResNorm.hpp
  BelosStatusTestLogResNorm.hpp
  BelosStatusTestMaxIters.hpp
  BelosStatusTestOutputFactory.hpp
  BelosStatusTestOutput.hpp
  BelosStatusTestResNorm.hpp
  BelosStatusTestResNormOutput.hpp
  BelosStatusTestUserOutput.hpp
  BelosStochasticCGIteration.hpp
  BelosStubTsqrAdapter.hpp
  BelosTFQMRIter.hpp
  BelosTFQMRSolMgr.hpp
  BelosPseudoBlockTFQMRIter.hpp
  BelosPseudoBlockTFQMRSolMgr.hpp
  BelosTypes.hpp
  Belos_Details_LinearSolver.hpp
  Belos_Details_LinearSolverFactory.hpp
  Belos_Details_registerLinearSolverFactory.hpp
  Belos_Details_registerSolverFactory.hpp
  BelosUtils.hpp
  )

# UTIL
APPEND_SET(SOURCES
  Belos_Details_EBelosSolverType.cpp
  BelosTypes.cpp
  BelosVersion.cpp
  Belos_Details_registerLinearSolverFactory.cpp
  BelosSolverFactory.cpp
  Belos_Details_registerSolverFactory.cpp
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
  BelosPseudoBlockTFQMRSolMgr.hpp
  BelosCGIter.hpp
  BelosCGSingleRedIter.hpp
  BelosBlockCGIter.hpp
  BelosMinresIter.hpp
  BelosBlockGmresIter.hpp
  BelosPseudoBlockCGIter.hpp
  BelosPseudoBlockGmresIter.hpp
  BelosPseudoBlockStochasticCGIter.hpp
  BelosTFQMRIter.hpp
  BelosPseudoBlockTFQMRIter.hpp
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
  BelosBlockGCRODRIter.hpp
  BelosBlockGCRODRSolMgr.hpp
  BelosInnerSolveResult.hpp
  BelosInnerSolver.hpp
  BelosProjectedLeastSquaresSolver.hpp
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
    )
ELSE ()
  TRIBITS_ADD_LIBRARY(
    belos
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )
ENDIF ()
