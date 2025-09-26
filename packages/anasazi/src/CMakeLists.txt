

#
# A) Package-specific configuration options
#

IF(Anasazi_ENABLE_ThyraCore)
  SET(HAVE_ANASAZI_THYRA 1)
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_TSQR)
IF(${PACKAGE_NAME}_ENABLE_TSQR)
  SET(HAVE_ANASAZI_TSQR 1)
ENDIF()

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(CURRENT_PACKAGE ANASAZI)
CONFIGURE_FILE("${Trilinos_SOURCE_DIR}/packages/Trilinos_DLLExportMacro.h.in"
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_DLLExportMacro.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_DLLExportMacro.h
  )

#
# Core Files
#

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

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
  AnasaziEigensolverDecl.hpp
  AnasaziEigensolver.hpp
  AnasaziFactory.hpp
  AnasaziGeneralizedDavidson.hpp
  AnasaziGeneralizedDavidsonSolMgr.hpp
  AnasaziGenOrthoManager.hpp
  AnasaziHelperTraits.hpp
  AnasaziICGSOrthoManager.hpp
  AnasaziIRTR.hpp
  AnasaziLOBPCG.hpp
  AnasaziLOBPCGSolMgr.hpp
  AnasaziMatOrthoManager.hpp
  AnasaziMinres.hpp
  AnasaziMultiVec.hpp
  AnasaziMultiVecTraits.hpp
  AnasaziOperator.hpp
  AnasaziOperatorTraits.hpp
  AnasaziOrthoManager.hpp
  AnasaziOutputManager.hpp
  AnasaziOutputStreamTraits.hpp
  AnasaziRTRBase.hpp
  AnasaziRTRSolMgr.hpp
  AnasaziSaddleContainer.hpp
  AnasaziSaddleOperator.hpp
  AnasaziSimpleLOBPCGSolMgr.hpp
  AnasaziSIRTR.hpp
  AnasaziSolverManager.hpp
  AnasaziSolverUtils.hpp
  AnasaziSortManager.hpp
  AnasaziStatusTestCombo.hpp
  AnasaziStatusTestDecl.hpp
  AnasaziStatusTest.hpp
  AnasaziStatusTestMaxIters.hpp
  AnasaziStatusTestOutput.hpp
  AnasaziStatusTestResNorm.hpp
  AnasaziStatusTestSpecTrans.hpp
  AnasaziStatusTestWithOrdering.hpp
  AnasaziStubTsqrAdapter.hpp
  AnasaziSVQBOrthoManager.hpp
  AnasaziTraceMinBase.hpp
  AnasaziTraceMinBaseSolMgr.hpp
  AnasaziTraceMinDavidson.hpp
  AnasaziTraceMinDavidsonSolMgr.hpp
  AnasaziTraceMin.hpp
  AnasaziTraceMinRitzOp.hpp
  AnasaziTraceMinSolMgr.hpp
  AnasaziTraceMinTypes.hpp
  AnasaziTypes.hpp
  AnasaziGlobalComm.hpp
  )

APPEND_SET(SOURCES
  AnasaziVersion.cpp
  )

#
# Util Files
#

APPEND_SET(HEADERS
  AnasaziMVOPTester.hpp
  )


#
# Comp Files
#

APPEND_SET(HEADERS
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
  AnasaziTraceMinBase.hpp
  AnasaziTraceMin.hpp
  AnasaziTraceMinDavidson.hpp
  AnasaziTraceMinBaseSolMgr.hpp
  AnasaziTraceMinSolMgr.hpp
  AnasaziTraceMinDavidsonSolMgr.hpp
  )

ASSERT_DEFINED(Anasazi_ENABLE_ThyraCore)
IF (Anasazi_ENABLE_ThyraCore)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../thyra/src)
  APPEND_SET(HEADERS
    ../thyra/src/AnasaziThyraAdapter.hpp
    ../thyra/src/AnasaziThyraDebugAdapter.hpp
    )
ENDIF()

ASSERT_DEFINED(Anasazi_ENABLE_Tpetra)
IF (Anasazi_ENABLE_Tpetra)
  TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../tpetra/src)
  # The AnasaziTpetraAdapter headers are already pulled in by
  # <./packages/anasazi/tpetra/src/CMakeLists.txt>.
ENDIF()

IF (TPL_ENABLE_MPI)

ENDIF()


######################################################################
## BEGIN TSQR setup
######################################################################

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_TSQR)
IF(${PACKAGE_NAME}_ENABLE_TSQR)
  #
  # Build the TSQR-based OrthoManager and MatOrthoManager classes.
  #
  APPEND_SET(HEADERS AnasaziTsqrOrthoManagerImpl.hpp)
  APPEND_SET(HEADERS AnasaziTsqrOrthoManager.hpp)
ENDIF() # ${PACKAGE_NAME}_ENABLE_TSQR

######################################################################
## END TSQR setup
######################################################################

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  anasazi
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEFINES -DANASAZI_LIB_EXPORTS_MODE
  )
