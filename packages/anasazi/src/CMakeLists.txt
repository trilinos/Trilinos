
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

# Enabling a define if the Thyra Epetra adapter is defined.
ASSERT_DEFINED(Thyra_ENABLE_ThyraEpetraAdapters)
IF(Thyra_ENABLE_ThyraEpetraAdapters)
  SET(HAVE_EPETRA_THYRA 1)
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_TSQR)
IF(${PACKAGE_NAME}_ENABLE_TSQR)
  SET(HAVE_ANASAZI_TSQR 1)
ENDIF()

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(CURRENT_PACKAGE ANASAZI)
CONFIGURE_FILE("${TRILINOS_HOME_DIR}/packages/Trilinos_DLLExportMacro.h.in"
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_DLLExportMacro.h)

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
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_DLLExportMacro.h
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
  )

ASSERT_DEFINED(Anasazi_ENABLE_Thyra)
IF (Anasazi_ENABLE_Thyra)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../thyra/src)
  APPEND_SET(HEADERS
    ../thyra/src/AnasaziThyraAdapter.hpp
    ../thyra/src/AnasaziThyraDebugAdapter.hpp
    )
ENDIF()

ASSERT_DEFINED(Anasazi_ENABLE_Tpetra)
IF (Anasazi_ENABLE_Tpetra)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../tpetra/src)
  APPEND_SET(HEADERS
    ../tpetra/src/AnasaziTpetraAdapter.hpp
    )
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

PACKAGE_ADD_LIBRARY(
  anasazi
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEFINES -DANASAZI_LIB_EXPORTS_MODE
  )
