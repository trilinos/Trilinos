
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

SET(EPETRA_MPI ${TPL_ENABLE_MPI})
IF(CMAKE_SIZEOF_VOID_P GREATER 4)
  SET(EPETRA_ADDRESS64BIT ON)
ENDIF()
SET(Epetra_HAVE_OMP ${Trilinos_ENABLE_OpenMP})
# This is the preferred macro to key on, but am keeping the old one for backward compatibility
SET(EPETRA_HAVE_OMP ${Trilinos_ENABLE_OpenMP})

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(CURRENT_PACKAGE EPETRA)
CONFIGURE_FILE("${TRILINOS_HOME_DIR}/packages/Trilinos_DLLExportMacro.h.in"
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_DLLExportMacro.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_DLLExportMacro.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  Epetra_BLAS.h
  Epetra_BLAS_wrappers.h
  Epetra_BasicDirectory.h
  Epetra_BasicRowMatrix.h
  Epetra_BlockMap.h
  Epetra_BlockMapData.h
  Epetra_C_wrappers.h
  Epetra_CombineMode.h
  Epetra_Comm.h
  Epetra_CompObject.h
  Epetra_ConfigDefs.h
  Epetra_CrsGraph.h
  Epetra_CrsGraphData.h
  Epetra_CrsMatrix.h
  Epetra_CrsSingletonFilter.h
  Epetra_Data.h
  Epetra_DataAccess.h
  Epetra_Directory.h
  Epetra_DistObject.h
  Epetra_Distributor.h
  Epetra_Export.h
  Epetra_FECrsGraph.h
  Epetra_FECrsMatrix.h
  Epetra_FEVbrMatrix.h
  Epetra_FEVector.h
#  Epetra_FastCrsMatrix.h
  Epetra_Flops.h
  Epetra_Fortran_wrappers.h
  Epetra_HashTable.h
  Epetra_Import.h
  Epetra_IntSerialDenseMatrix.h
  Epetra_IntSerialDenseVector.h
  Epetra_IntVector.h
  Epetra_InvOperator.h
  Epetra_JadMatrix.h
  Epetra_LAPACK.h
  Epetra_LAPACK_wrappers.h
  Epetra_LinearProblem.h
#  Epetra_LinearProblemRedistor.h
  Epetra_LocalMap.h
  Epetra_Map.h
  Epetra_MapColoring.h
  Epetra_MultiVector.h
  Epetra_Object.h
  Epetra_OffsetIndex.h
  Epetra_Operator.h
  Epetra_RowMatrix.h
  Epetra_RowMatrixTransposer.h
  Epetra_SerialComm.h
  Epetra_SerialCommData.h
  Epetra_SerialDenseMatrix.h
  Epetra_SerialDenseOperator.h
  Epetra_SerialDenseSVD.h
  Epetra_SerialDenseSolver.h
  Epetra_SerialDenseVector.h
  Epetra_SerialDistributor.h
  Epetra_SerialSpdDenseSolver.h
  Epetra_SerialSymDenseMatrix.h
  Epetra_SrcDistObject.h
  Epetra_Time.h
  Epetra_Util.h
  Epetra_VbrMatrix.h
  Epetra_VbrRowMatrix.h
  Epetra_Vector.h
  Epetra_Version.h
  )

# Set sources
SET(SOURCES ${SOURCES}
  Epetra_BLAS.cpp
  Epetra_BasicDirectory.cpp
  Epetra_BasicRowMatrix.cpp
  Epetra_BlockMap.cpp
  Epetra_BlockMapData.cpp
  Epetra_C_wrappers.cpp
  Epetra_CompObject.cpp
  Epetra_CrsGraph.cpp
  Epetra_CrsGraphData.cpp
  Epetra_CrsMatrix.cpp
  Epetra_CrsSingletonFilter.cpp
  Epetra_Data.cpp
  Epetra_DistObject.cpp
  Epetra_Export.cpp
  Epetra_FECrsGraph.cpp
  Epetra_FECrsMatrix.cpp
  Epetra_FEVbrMatrix.cpp
  Epetra_FEVector.cpp
#  Epetra_FastCrsMatrix.cpp
  Epetra_Flops.cpp
  Epetra_Fortran_wrappers.cpp
  Epetra_Import.cpp
  Epetra_IntSerialDenseMatrix.cpp
  Epetra_IntSerialDenseVector.cpp
  Epetra_IntVector.cpp
#  Epetra_InvOperator.cpp
  Epetra_JadMatrix.cpp
  Epetra_LAPACK.cpp
  Epetra_LinearProblem.cpp
#  Epetra_LinearProblemRedistor.cpp
  Epetra_LocalMap.cpp
  Epetra_Map.cpp
  Epetra_MapColoring.cpp
  Epetra_MultiVector.cpp
  Epetra_Object.cpp
  Epetra_OffsetIndex.cpp
  Epetra_RowMatrixTransposer.cpp
  Epetra_SerialComm.cpp
  Epetra_SerialCommData.cpp
  Epetra_SerialDenseMatrix.cpp
  Epetra_SerialDenseSVD.cpp
  Epetra_SerialDenseSolver.cpp
  Epetra_SerialDenseVector.cpp
  Epetra_SerialDistributor.cpp
  Epetra_SerialSpdDenseSolver.cpp
  Epetra_SerialSymDenseMatrix.cpp
  Epetra_Time.cpp
  Epetra_Util.cpp
  Epetra_VbrMatrix.cpp
  Epetra_Vector.cpp
  )

IF(Epetra_ENABLE_Oski)
  APPEND_SET(HEADERS
    Epetra_OskiUtils.h
    Epetra_OskiMultiVector.h
    Epetra_OskiVector.h
    Epetra_OskiError.h
    Epetra_OskiPermutation.h
    Epetra_OskiMatrix.h
    )

  APPEND_SET(SOURCES
    Epetra_OskiUtils.cpp
    Epetra_OskiMultiVector.cpp
    Epetra_OskiVector.cpp
    Epetra_OskiError.cpp
    Epetra_OskiPermutation.cpp
    Epetra_OskiMatrix.cpp
    )
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Fortran)
  SET(SOURCES ${SOURCES}	
    Epetra_dcrsmm.f
    Epetra_dcrsmv.f
    Epetra_dcrssm.f
    Epetra_dcrssv.f
    )	
ENDIF()

# Set addition headers and source based on package-specific 
# configuration options

IF(TPL_ENABLE_MPI)

  SET(HEADERS ${HEADERS}
    Epetra_MpiComm.h
    Epetra_MpiCommData.h
    Epetra_MpiDistributor.h
    )

  SET(SOURCES ${SOURCES}
    Epetra_MpiComm.cpp
    Epetra_MpiCommData.cpp
    Epetra_MpiDistributor.cpp
    )

ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  epetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEFINES -DEPETRA_LIB_EXPORTS_MODE
  )
