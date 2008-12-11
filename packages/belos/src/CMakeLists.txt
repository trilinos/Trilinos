
INCLUDE(TrilinosPackageLibraryHelpers)

#
# A) Package-specific configuration options
#

TRILINOS_PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

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
SET(HEADERS ${HEADERS}
  BelosBlockCGSolMgr.hpp
  BelosBlockCGIter.hpp
  BelosBlockGmresIter.hpp
  BelosBlockFGmresIter.hpp
  BelosBlockGmresSolMgr.hpp
  BelosCGIter.hpp
  BelosCGIteration.hpp
  BelosConfigDefs.hpp
  BelosDGKSOrthoManager.hpp
  BelosEpetraAdapter.hpp
  BelosEpetraOperator.h
  BelosPCPGIter.hpp
  BelosPCPGSolMgr.hpp
  BelosGCRODRIter.hpp
  BelosGCRODRSolMgr.hpp
  BelosGmresIteration.hpp
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
  BelosStatusTestResNorm.hpp
  BelosStatusTestGenResNorm.hpp
  BelosStatusTestImpResNorm.hpp
  BelosTypes.hpp
  ) 

# CORE
SET(SOURCES ${SOURCES}
  BelosEpetraAdapter.cpp
  BelosEpetraOperator.cpp
  )

# UTIL
SET(SOURCES ${SOURCES}
  BelosVersion.cpp
  )

# COMP_H
SET(HEADERS ${HEADERS}
  BelosBlockCGSolMgr.hpp
  BelosBlockGmresSolMgr.hpp
  BelosPseudoBlockGmresSolMgr.hpp
  BelosCGIter.hpp
  BelosBlockCGIter.hpp
  BelosBlockGmresIter.hpp
  BelosPseudoBlockGmresIter.hpp
  )

#
# C) Define the targets for package's library(s)
#

TRILINOS_PACKAGE_ADD_LIBRARY(
  belos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

#
# D) Export the dependancy variables of this package for other packages to use
#

TRILINOS_PACKAGE_EXPORT_DEPENDENCY_VARIABLES()
