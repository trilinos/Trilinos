
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

# Enabling a define if the Thyra Epetra adapter is defined.
ASSERT_DEFINED(Thyra_ENABLE_Epetra)
IF(Thyra_ENABLE_Epetra)
  SET(HAVE_EPETRA_THYRA 1)
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
## mfh 28 Jun 2010
## BEGIN TSQR setup
######################################################################

# Enable Fortran (used by TSQR, for one file)
enable_language (Fortran)
# Generate macros for Fortran name mangling, for use by TSQR
include(FortranCInterface)
FortranCInterface_HEADER("TSQR/Tsqr_FortranCInterface.hpp")
# Include tests for LAPACK functionality required by TSQR.
include ("TSQR/TsqrTestLapack.cmake")
# Generate TSQR's LAPACK wrappers.
include ("TSQR/TsqrLapackWrappers.cmake")

# Include headers and source files for Trilinos -> TSQR adaptors.
# Some of these depend on Epetra or Tpetra, so check accordingly if
# Anasazi_ENABLE_{Epetra,Tpetra} is defined.
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/TSQR)
APPEND_SET(HEADERS 
  TSQR/TsqrAdaptor.hpp
  TSQR/TsqrCommFactory.hpp
  TSQR/TsqrFactory.hpp
  TSQR/TsqrFactory_SequentialTsqr.hpp
  TSQR/TsqrFactory_TbbTsqr.hpp
  TSQR/TsqrRandomizer.hpp
  TSQR/TsqrTypeAdaptor.hpp
  )
APPEND_SET(HEADERS
  TSQR/Tsqr_ApplyType.hpp
  TSQR/Tsqr_Blas.hpp
  TSQR/Tsqr_CacheBlocker.hpp
  TSQR/Tsqr_CacheBlockingStrategy.hpp
  TSQR/Tsqr_Combine.hpp
  TSQR/Tsqr_Config.hpp
  TSQR/Tsqr_DistTsqrHelper.hpp
  TSQR/Tsqr_DistTsqr.hpp
  TSQR/Tsqr_FortranCInterface.hpp
  TSQR/Tsqr_generateStack.hpp
  TSQR/Tsqr_GlobalVerify.hpp
  TSQR/Tsqr.hpp
  TSQR/Tsqr_Lapack.hpp
  TSQR/Tsqr_LocalVerify.hpp
  TSQR/Tsqr_Matrix.hpp
  TSQR/Tsqr_MatView.hpp
  TSQR/Tsqr_MessengerBase.hpp
  TSQR/Tsqr_Mgs.hpp
  TSQR/Tsqr_MgsTest.hpp
  TSQR/Tsqr_nodeTestProblem.hpp
  TSQR/Tsqr_NodeTsqr.hpp
  TSQR/Tsqr_ParTest.hpp
  TSQR/Tsqr_printGlobalMatrix.hpp
  TSQR/Tsqr_Random_GlobalMatrix.hpp
  TSQR/Tsqr_Random_MatrixGenerator.hpp
  TSQR/Tsqr_Random_NormalGenerator.hpp
  TSQR/Tsqr_RMessenger.hpp
  TSQR/Tsqr_ScalarTraits.hpp
  TSQR/Tsqr_SeqTest.hpp
  TSQR/Tsqr_SequentialCholeskyQR.hpp
  TSQR/Tsqr_SequentialTsqr.hpp
  TSQR/Tsqr_TbbTest.hpp
  TSQR/Tsqr_TestSetup.hpp
  TSQR/Tsqr_TsqrTest.hpp
  TSQR/Tsqr_Util.hpp
  TSQR/Tsqr_verifyTimerConcept.hpp
  )
APPEND_SET(SOURCES
  TSQR/Tsqr_ApplyType.cpp 
  TSQR/Tsqr_Blas.cpp 
  TSQR/Tsqr_DLapack.cpp
  TSQR/Tsqr_SLapack.cpp
  TSQR/Tsqr_ZLapack.cpp
  TSQR/Tsqr_CLapack.cpp
  TSQR/Tsqr_Combine.cpp
  TSQR/Tsqr_HouseholderReflector.f90 # generated by TsqrLapackWrappers.cmake (in this directory)
  TSQR/Tsqr_LapackQrRoutines.f90
  TSQR/Tsqr_Combine.f90
  )
APPEND_SET(HEADERS
  TSQR/TsqrFactory.hpp
  TSQR/TsqrFactory_SequentialTsqr.hpp
  )
#IF (Anasazi_ENABLE_Epetra)
#  APPEND_SET(HEADERS
#    TsqrAdaptor_Epetra_MultiVector.hpp
#    TsqrCommFactory_Epetra.hpp
#    TsqrRandomizer_Epetra_MultiVector.hpp
#    Tsqr_EpetraMessenger.hpp # FIXME (mfh 15 Jul 2010) This file doesn't exist yet.
#    TsqrTypeAdaptor_Epetra_MultiVector.hpp
#    )
#ENDIF()

ASSERT_DEFINED(Anasazi_ENABLE_Tpetra)
IF (Anasazi_ENABLE_Tpetra)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../tpetra/src)
  APPEND_SET(HEADERS
    ../tpetra/src/TsqrAdaptor_Tpetra_MultiVector.hpp
    ../tpetra/src/TsqrCommFactory_Tpetra.hpp
    ../tpetra/src/TsqrRandomizer_Tpetra_MultiVector.hpp
    ../tpetra/src/Tsqr_TpetraMessenger.hpp
    ../tpetra/src/TsqrTypeAdaptor_Tpetra_MultiVector_SerialNode.hpp
    ../tpetra/src/TsqrTypeAdaptor_Tpetra_MultiVector_TBBNode.hpp
    )
  # FIXME (mfh 29 Jun 2010) Should test for whether we are linking
  # with TBB, and only append these headers if we are.
  APPEND_SET(HEADERS
    TSQR/TbbTsqr_ApplyTask.hpp
    TSQR/TbbTsqr_CacheBlockTask.hpp
    TSQR/TbbTsqr_ExplicitQTask.hpp
    TSQR/TbbTsqr_FactorTask.hpp
    TSQR/TbbTsqr_FillWithZerosTask.hpp
    TSQR/TbbTsqr.hpp
    TSQR/TbbTsqr_Partitioner.hpp
    TSQR/TbbTsqr_TbbMgs.hpp
    TSQR/TbbTsqr_TbbParallelTsqr.hpp
    TSQR/TbbTsqr_UnCacheBlockTask.hpp
    TSQR/TsqrFactory_TbbTsqr.hpp
    )
ENDIF()

######################################################################
## mfh 28 Jun 2010
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
