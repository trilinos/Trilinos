
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
  BelosBlockCGSolMgr.hpp
  BelosBlockCGIter.hpp
  BelosBlockGmresIter.hpp
  BelosBlockFGmresIter.hpp
  BelosBlockGmresSolMgr.hpp
  BelosCGIter.hpp
  BelosCGIteration.hpp
  BelosConfigDefs.hpp
  BelosDGKSOrthoManager.hpp
  BelosPCPGIter.hpp
  BelosPCPGSolMgr.hpp
  BelosGCRODRIter.hpp
  BelosGCRODRSolMgr.hpp
  BelosGmresIteration.hpp
  BelosGmresPolyOp.hpp
  BelosGmresPolySolMgr.hpp
  BelosICGSOrthoManager.hpp
  BelosIMGSOrthoManager.hpp
  BelosIteration.hpp
  BelosLinearProblem.hpp
  BelosMatOrthoManager.hpp
  BelosMultiVec.hpp
  BelosMultiVecTraits.hpp
  BelosMVOPTester.hpp
  BelosOperator.hpp
  BelosOperatorTraits.hpp
  BelosOrthoManager.hpp
  BelosOutputManager.hpp
  BelosPseudoBlockCGIter.hpp
  BelosPseudoBlockCGSolMgr.hpp
  BelosPseudoBlockGmresIter.hpp
  BelosPseudoBlockGmresSolMgr.hpp
  BelosRCGIter.hpp
  BelosRCGSolMgr.hpp
  BelosSolverManager.hpp
  BelosStatusTest.hpp
  BelosStatusTestCombo.hpp
  BelosStatusTestMaxIters.hpp
  BelosStatusTestOutput.hpp
  BelosStatusTestOutputFactory.hpp
  BelosStatusTestGeneralOutput.hpp
  BelosStatusTestResNormOutput.hpp
  BelosStatusTestResNorm.hpp
  BelosStatusTestGenResNorm.hpp
  BelosStatusTestImpResNorm.hpp
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
  BelosBlockGmresSolMgr.hpp
  BelosPseudoBlockCGSolMgr.hpp
  BelosPseudoBlockGmresSolMgr.hpp
  BelosTFQMRSolMgr.hpp
  BelosCGIter.hpp
  BelosBlockCGIter.hpp
  BelosBlockGmresIter.hpp
  BelosPseudoBlockCGIter.hpp
  BelosPseudoBlockGmresIter.hpp
  BelosTFQMRIter.hpp
  )

ASSERT_DEFINED(Belos_ENABLE_Tpetra)
IF (Belos_ENABLE_Tpetra)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../tpetra/src)
  APPEND_SET(HEADERS
    ../tpetra/src/BelosTpetraAdapter.hpp
    )
ENDIF()


#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  belos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
