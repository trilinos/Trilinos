
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
FortranCInterface_HEADER("Tsqr_FortranCInterface.hpp")
# Include tests for LAPACK functionality required by TSQR.
include ("TsqrTestLapack.cmake")

# Include headers and source files for Trilinos -> TSQR adaptors.
# Some of these depend on Epetra or Tpetra, so check accordingly if
# Anasazi_ENABLE_{Epetra,Tpetra} is defined.
APPEND_SET(HEADERS 
  TsqrAdaptor.hpp
  TsqrFactory.hpp
  TsqrFactory_SequentialTsqr.hpp
  TsqrFactory_TbbTsqr.hpp
  TsqrTrilinosMessenger.hpp
  TsqrTypeAdaptor.hpp
  )
APPEND_SET(HEADERS
  Tsqr_ApplyType.hpp
  Tsqr_Blas.hpp
  Tsqr_CacheBlocker.hpp
  Tsqr_CacheBlockingStrategy.hpp
  Tsqr_Combine.hpp
  Tsqr_DistTsqrHelper.hpp
  Tsqr_DistTSqr.hpp
  Tsqr_GlobalVerify.hpp
  Tsqr.hpp
  Tsqr_Lapack.hpp
  Tsqr_LocalVerify.hpp
  Tsqr_Matrix.hpp
  Tsqr_MatView.hpp
  Tsqr_Mgs.hpp
  Tsqr_NodeTsqr.hpp
  Tsqr_Random_MatrixGenerator.hpp
  Tsqr_ScalarTraits.hpp
  Tsqr_SequentialCholeskyQR.hpp
  Tsqr_SequentialTsqr.hpp
  Tsqr_Util.hpp
  )
APPEND_SET(SOURCES
  Tsqr_ApplyType.cpp 
  Tsqr_Blas.cpp 
  Tsqr_DLapack.cpp
  Tsqr_SLapack.cpp
  Tsqr_ZLapack.cpp
  Tsqr_CLapack.cpp
  Tsqr_Combine.cpp
  Tsqr_LapackQrRoutines.f90
  Tsqr_Combine.f90
  )
IF (Anasazi_ENABLE_Epetra)
  APPEND_SET(HEADERS
    TsqrAdaptor_Epetra_MultiVector.hpp
    TsqrTypeAdaptor_Epetra_MultiVector.hpp
    )
ENDIF()
IF (Anasazi_ENABLE_Tpetra)
  APPEND_SET(HEADERS
    TsqrAdaptor_Tpetra_MultiVector_SerialNode.hpp
    TsqrAdaptor_Tpetra_MultiVector_TBBNode.hpp
    TsqrTypeAdaptor_Tpetra_MultiVector_SerialNode.hpp
    TsqrTypeAdaptor_Tpetra_MultiVector_TBBNode.hpp
    )
  # FIXME (mfh 29 Jun 2010) Should test for whether we are linking
  # with TBB, and only append these headers if we are.
  APPEND_SET(HEADERS
    TbbTsqr_ApplyTask.hpp
    TbbTsqr_CacheBlockTask.hpp
    TbbTsqr_ExplicitQTask.hpp
    TbbTsqr_FactorTask.hpp
    TbbTsqr_FillWithZerosTask.hpp
    TbbTsqr.hpp
    TbbTsqr_Partitioner.hpp
    TbbTsqr_TbbMgs.hpp
    TbbTsqr_TbbParallelTsqr.hpp
    TbbTsqr_UnCacheBlockTask.hpp
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
