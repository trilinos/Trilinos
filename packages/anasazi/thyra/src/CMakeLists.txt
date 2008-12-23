
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

#
# Core Files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  AnasaziBasicEigenproblem.hpp
  AnasaziBasicOrthoManager.hpp
  AnasaziBasicOutputManager.hpp
  AnasaziBasicSort.hpp
  AnasaziBlockDavidson.hpp
  AnasaziBlockDavidsonSolMgr.hpp
  AnasaziBlockKrylovSchur.hpp
  AnasaziBlockKrylovSchurSolMgr.hpp
  AnasaziConfigDefs.hpp
  AnasaziEigenproblem.hpp
  AnasaziEigensolver.hpp
  AnasaziEigensolverDecl.hpp
  AnasaziLOBPCG.hpp
  AnasaziLOBPCGSolMgr.hpp
  AnasaziGenOrthoManager.hpp
  AnasaziHelperTraits.hpp
  AnasaziICGSOrthoManager.hpp
  AnasaziIRTR.hpp
  AnasaziMatOrthoManager.hpp
  AnasaziMultiVec.hpp
  AnasaziMultiVecTraits.hpp
  AnasaziOperator.hpp
  AnasaziOperatorTraits.hpp
  AnasaziOrthoManager.hpp
  AnasaziOutputManager.hpp
  AnasaziRTRBase.hpp
  AnasaziRTRSolMgr.hpp
  AnasaziSimpleLOBPCGSolMgr.hpp
  AnasaziSIRTR.hpp
  AnasaziSolverManager.hpp
  AnasaziSortManager.hpp
  AnasaziSolverUtils.hpp
  AnasaziStatusTest.hpp
  AnasaziStatusTestDecl.hpp
  AnasaziStatusTestCombo.hpp
  AnasaziStatusTestMaxIters.hpp
  AnasaziStatusTestWithOrdering.hpp
  AnasaziStatusTestOutput.hpp
  AnasaziStatusTestResNorm.hpp
  AnasaziSVQBOrthoManager.hpp
  AnasaziTypes.hpp
  )

APPEND_SET(SOURCES ${SOURCES}
  AnasaziVersion.cpp
  )

#
# Util Files
#

APPEND_SET(HEADERS ${HEADERS}
  AnasaziMVOPTester.hpp
  )


#
# Comp Files
#

APPEND_SET(HEADERS ${HEADERS}
  AnasaziBlockDavidson.hpp
  AnasaziBlockKrylovSchur.hpp
  AnasaziLOBPCG.hpp
  AnasaziRTRBase.hpp
  AnasaziIRTR.hpp
  AnasaziSIRTR.hpp
  AnasaziBlockDavidsonSolMgr.hpp
  AnasaziBlockKrylovSchurSolMgr.hpp
  AnasaziLOBPCGSolMgr.hpp
  AnasaziSimpleLOBPCGSolMgr.hpp
  AnasaziRTRSolMgr.hpp
  )

IF (TPL_ENABLE_MPI)

ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  anasazi
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
