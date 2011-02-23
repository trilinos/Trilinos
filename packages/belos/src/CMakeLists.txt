
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# CORE
APPEND_SET(HEADERS
  BelosAkxFactory.hpp
  BelosAkx.hpp
  BelosBlockCGIter.hpp
  BelosBlockCGSolMgr.hpp
  BelosBlockFGmresIter.hpp
  BelosBlockGmresIter.hpp
  BelosBlockGmresSolMgr.hpp
  BelosCaGmres.hpp  # Subclass of GmresBase
  BelosCGIteration.hpp
  BelosCGIter.hpp
  BelosConfigDefs.hpp
  BelosDGKSOrthoManager.hpp
  BelosGCRODRIter.hpp
  BelosGCRODRSolMgr.hpp
  BelosGmresBaseFactory.hpp
  BelosGmresBase.hpp
  BelosGmresBaseIteration.hpp
  BelosGmresIteration.hpp
  BelosGmresPolyOp.hpp
  BelosGmresPolySolMgr.hpp
  BelosGmresSolMgr.hpp
  BelosICGSOrthoManager.hpp
  BelosIMGSOrthoManager.hpp
  BelosIteration.hpp
  BelosLinearProblem.hpp
  BelosMatOrthoManager.hpp
  BelosMinresIteration.hpp
  BelosMinresIter.hpp
  BelosMinresSolMgr.hpp
  BelosMonomialOpAkx.hpp
  BelosMultiVec.hpp
  BelosMultiVecTraits.hpp
  BelosMVOPTester.hpp
  BelosNewtonOpAkx.hpp
  BelosOpAkx.hpp
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
  BelosRCGIter.hpp
  BelosRCGSolMgr.hpp
  BelosSimpleOrthoManager.hpp
  BelosSolverManager.hpp
  BelosStandardGmres.hpp
  BelosStatusTestCombo.hpp
  BelosStatusTestGeneralOutput.hpp
  BelosStatusTestGenResNorm.hpp
  BelosStatusTest.hpp
  BelosStatusTestImpResNorm.hpp
  BelosStatusTestMaxIters.hpp
  BelosStatusTestOutputFactory.hpp
  BelosStatusTestOutput.hpp
  BelosStatusTestResNorm.hpp
  BelosStatusTestResNormOutput.hpp
  BelosTFQMRIter.hpp
  BelosTFQMRSolMgr.hpp
  BelosTypes.hpp
  ) 

# CORE
APPEND_SET(SOURCES
  )

# UTIL
APPEND_SET(SOURCES
  BelosVersion.cpp
  )

# COMP_H
APPEND_SET(HEADERS
  BelosBlockCGSolMgr.hpp
  BelosMinresSolMgr.hpp
  BelosBlockGmresSolMgr.hpp
  BelosPseudoBlockCGSolMgr.hpp
  BelosPseudoBlockGmresSolMgr.hpp
  BelosTFQMRSolMgr.hpp
  BelosCGIter.hpp
  BelosBlockCGIter.hpp
  BelosMinresIter.hpp
  BelosBlockGmresIter.hpp
  BelosPseudoBlockCGIter.hpp
  BelosPseudoBlockGmresIter.hpp
  BelosTFQMRIter.hpp
  )

######################################################################
## BEGIN TSQR setup
######################################################################

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_TSQR)
IF(${PACKAGE_NAME}_ENABLE_TSQR)
  #
  # Build the TSQR-based OrthoManager and MatOrthoManager classes.
  #
  APPEND_SET(HEADERS BelosTsqrOrthoManagerImpl.hpp BelosTsqrOrthoManager.hpp)
ENDIF() # ${PACKAGE_NAME}_ENABLE_TSQR

######################################################################
## END TSQR setup
######################################################################

#
# C) Define the targets for package's library(s)
#

IF (${PACKAGE_NAME}_ENABLE_TSQR)
  PACKAGE_ADD_LIBRARY(
    belos
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS kokkos kokkosnodetsqr kokkosdisttsqr
    )
ELSE ()
  PACKAGE_ADD_LIBRARY(
    belos
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )
ENDIF ()
