
INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_PETSC)
SET(HAVE_PETSC ${${PACKAGE_NAME}_ENABLE_PETSC})
SET(HAVE_LIBPARMETIS ${Zoltan_ENABLE_ParMETIS})

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

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

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  EpetraExt_ConfigDefs.h
  EpetraExt_MMHelpers.h
  EpetraExt_MatrixMatrix.h
  EpetraExt_ProductOperator.h
  EpetraExt_TimedEpetraOperator.hpp
  EpetraExt_Version.h
  )

SET(SOURCES ${SOURCES}
  EpetraExt_ProductOperator.cpp
  EpetraExt_MMHelpers.cpp
  EpetraExt_MatrixMatrix.cpp
  EpetraExt_MatrixMatrix_mult_A_B.cpp
  EpetraExt_TimedEpetraOperator.cpp
  )

#
# src/transform
#

IF (${PACKAGE_NAME}_BUILD_TRANSFORM)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/transform)

  SET(HEADERS ${HEADERS}
    transform/EpetraExt_SubCopy_CrsMatrix.h
    transform/EpetraExt_Transform.h
    transform/EpetraExt_Transform_Composite.h
    transform/EpetraExt_Permutation.h
    transform/EpetraExt_Permutation_impl.h
    transform/EpetraExt_LPTrans_From_GraphTrans.h
    transform/EpetraExt_LPTrans_From_MatrixTrans.h
    transform/EpetraExt_View_CrsGraph.h
    transform/EpetraExt_View_CrsMatrix.h
    transform/EpetraExt_View_MultiVector.h
    transform/EpetraExt_Transpose_CrsGraph.h
    transform/EpetraExt_Transpose_RowMatrix.h
    transform/EpetraExt_Scale_LinearProblem.h
    transform/EpetraExt_CrsSingletonFilter_LinearProblem.h
    transform/EpetraExt_Overlap_CrsGraph.h
    transform/EpetraExt_Reindex_CrsMatrix.h
    transform/EpetraExt_Reindex_MultiVector.h
    transform/EpetraExt_Reindex_LinearProblem.h
    transform/EpetraExt_Reindex_LinearProblem2.h
    transform/EpetraExt_SolverMap_CrsMatrix.h
    transform/EpetraExt_SolverMap_LinearProblem.h
    transform/EpetraExt_BlockAdjacencyGraph.h
    )

  SET(SOURCES ${SOURCES}
    transform/EpetraExt_SubCopy_CrsMatrix.cpp
    transform/EpetraExt_Permutation.cpp
    transform/EpetraExt_LPTrans_From_GraphTrans.cpp
    transform/EpetraExt_LPTrans_From_MatrixTrans.cpp
    transform/EpetraExt_View_CrsGraph.cpp
    transform/EpetraExt_View_CrsMatrix.cpp
    transform/EpetraExt_View_MultiVector.cpp
    transform/EpetraExt_Transpose_CrsGraph.cpp
    transform/EpetraExt_Transpose_RowMatrix.cpp
    transform/EpetraExt_Scale_LinearProblem.cpp
    transform/EpetraExt_CrsSingletonFilter_LinearProblem.cpp
    transform/EpetraExt_Overlap_CrsGraph.cpp
    transform/EpetraExt_Reindex_CrsMatrix.cpp
    transform/EpetraExt_Reindex_MultiVector.cpp
    transform/EpetraExt_Reindex_LinearProblem.cpp
    transform/EpetraExt_Reindex_LinearProblem2.cpp
    transform/EpetraExt_SolverMap_CrsMatrix.cpp
    transform/EpetraExt_SolverMap_LinearProblem.cpp
    transform/EpetraExt_BlockAdjacencyGraph.cpp
    )

  IF (${PACKAGE_NAME}_BUILD_GRAPH_REORDERINGS AND ${PACKAGE_NAME}_BUILD_EXPERIMENTAL)

    SET(HEADERS ${HEADERS}
      transform/EpetraExt_SymmRCM_CrsGraph.h
      transform/EpetraExt_AMD_CrsGraph.h
      )

    SET(SOURCES ${SOURCES}
      transform/EpetraExt_SymmRCM_CrsGraph.cpp
      transform/EpetraExt_AMD_CrsGraph.cpp
      )

  ENDIF()

ENDIF()

#
# src/btf
#

IF (${PACKAGE_NAME}_BUILD_BTF)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/btf)

  SET(HEADERS ${HEADERS}
    btf/EpetraExt_BTF_CrsGraph.h
    btf/EpetraExt_BTF_CrsMatrix.h
    btf/EpetraExt_BTF_LinearProblem.h
    btf/EpetraExt_BlockJacobi_LinearProblem.h
    )

  SET(SOURCES ${SOURCES}
    btf/EpetraExt_BTF_CrsGraph.cpp
    btf/EpetraExt_BTF_CrsMatrix.cpp
    btf/EpetraExt_BTF_LinearProblem.cpp
    btf/EpetraExt_BlockJacobi_LinearProblem.cpp
    btf/pothen/concmp.f
    btf/pothen/genbtf.f
    btf/pothen/corsum.f
    btf/pothen/finsum.f
    btf/pothen/fnrsum.f
    btf/pothen/izero.f
    btf/pothen/maxmatch.f
    btf/pothen/rectblk.f
    btf/pothen/mmc13e.f
    btf/pothen/prtivs.f
    btf/pothen/strchk.f
    btf/pothen/mattrans.f
    )

ENDIF()

#
# src/inout
#

IF (${PACKAGE_NAME}_BUILD_INOUT)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/inout)

  SET(HEADERS ${HEADERS}
    inout/EpetraExt_mmio.h
    inout/EpetraExt_OperatorOut.h
    inout/EpetraExt_RowMatrixOut.h
    inout/EpetraExt_MultiVectorOut.h
    inout/EpetraExt_VectorOut.h
    inout/EpetraExt_BlockMapOut.h
    inout/EpetraExt_BlockMapIn.h
    inout/EpetraExt_CrsMatrixIn.h
    inout/EpetraExt_MultiVectorIn.h
    inout/EpetraExt_VectorIn.h
    inout/EpetraExt_Exception.h
    inout/EpetraExt_Utils.h
    inout/EpetraExt_XMLReader.h
    inout/EpetraExt_XMLWriter.h
    inout/EpetraExt_HDF5.h
    inout/EpetraExt_HDF5_Handle.h
    inout/EpetraExt_DistArray.h
    )

  SET(SOURCES ${SOURCES}
    inout/EpetraExt_mmio.cpp
    inout/EpetraExt_OperatorOut.cpp
    inout/EpetraExt_RowMatrixOut.cpp
    inout/EpetraExt_MultiVectorOut.cpp
    inout/EpetraExt_VectorOut.cpp
    inout/EpetraExt_BlockMapOut.cpp
    inout/EpetraExt_BlockMapIn.cpp
    inout/EpetraExt_CrsMatrixIn.cpp
    inout/EpetraExt_MultiVectorIn.cpp
    inout/EpetraExt_VectorIn.cpp
    inout/EpetraExt_Utils.cpp
    inout/EpetraExt_XMLReader.cpp
    inout/EpetraExt_XMLWriter.cpp
    inout/EpetraExt_HDF5.cpp
    inout/EpetraExt_HDF5_DistObject.cpp
    )

  IF (${PACKAGE_NAME}_ENABLE_Triutils)

    SET(HEADERS ${HEADERS}
      inout/EpetraExt_readEpetraLinearSystem.h
      )

    SET(SOURCES ${SOURCES}
      inout/EpetraExt_readEpetraLinearSystem.cpp
      )

  ENDIF()

ENDIF()

#
# src/zoltan
#

IF (${PACKAGE_NAME}_BUILD_Zoltan)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/zoltan)

  SET(HEADERS ${HEADERS}
    zoltan/EpetraExt_ZoltanQuery.h
    zoltan/Zoltan_LoadBalance.h
    zoltan/Zoltan_MigrationContainer.h
    zoltan/Zoltan_MigrationFunctions.h
    zoltan/Zoltan_MigrationObject.h
    zoltan/Zoltan_QueryContainer.h
    zoltan/Zoltan_QueryFunctions.h
    zoltan/Zoltan_QueryObject.h
    )

  SET(SOURCES ${SOURCES}
    zoltan/EpetraExt_ZoltanQuery.cpp
    zoltan/Zoltan_LoadBalance.cpp
    zoltan/Zoltan_MigrationContainer.cpp
    zoltan/Zoltan_MigrationFunctions.cpp
    zoltan/Zoltan_MigrationObject.cpp
    zoltan/Zoltan_QueryContainer.cpp
    zoltan/Zoltan_QueryFunctions.cpp
    zoltan/Zoltan_QueryObject.cpp
    )

  IF(TPL_ENABLE_MPI)
    SET(HEADERS ${HEADERS}
      zoltan/EpetraExt_Zoltan_CrsGraph.h
      zoltan/EpetraExt_ZoltanOrder_CrsGraph.h
      )
    SET(SOURCES ${SOURCES}
      zoltan/EpetraExt_Zoltan_CrsGraph.cpp
      zoltan/EpetraExt_ZoltanOrder_CrsGraph.cpp
      )
  ENDIF()
ENDIF()

#
# src/coloring
#

IF (${PACKAGE_NAME}_BUILD_COLORING)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/coloring)

  SET(HEADERS ${HEADERS}
    coloring/EpetraExt_MapColoring.h
    coloring/EpetraExt_MapColoringIndex.h
    )

  SET(SOURCES ${SOURCES}
    coloring/EpetraExt_MapColoring.cpp
    coloring/EpetraExt_MapColoringIndex.cpp
    )

ENDIF()

#
# src/matlab
#

IF (${PACKAGE_NAME}_BUILD_MATLAB)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/matlab)

  SET(HEADERS ${HEADERS}
    matlab/EpetraExt_MatlabEngine.h
    matlab/EpetraExt_PutMultiVector.h
    matlab/EpetraExt_PutRowMatrix.h
    matlab/EpetraExt_PutBlockMap.h
    )

  SET(SOURCES ${SOURCES}
    matlab/EpetraExt_MatlabEngine.cpp
    matlab/EpetraExt_PutMultiVector.cpp
    matlab/EpetraExt_PutRowMatrix.cpp
    matlab/EpetraExt_PutBlockMap.cpp
    )

ENDIF()

#
# src/model_evaluator
#

IF (${PACKAGE_NAME}_BUILD_MODEL_EVALUATOR)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/model_evaluator)

  SET(HEADERS ${HEADERS}
    model_evaluator/EpetraExt_ModelEvaluator.h
    model_evaluator/EpetraExt_ModelEvaluatorScalingTools.h
    model_evaluator/EpetraExt_PolynomialVectorTraits.h
    )

  SET(SOURCES ${SOURCES}
    model_evaluator/EpetraExt_ModelEvaluator.cpp
    model_evaluator/EpetraExt_ModelEvaluatorScalingTools.cpp
    )

ENDIF()

#
# src/block
#

IF (${PACKAGE_NAME}_BUILD_BLOCK)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/block)

  SET(HEADERS ${HEADERS}
    block/EpetraExt_BlockCrsMatrix.h
    block/EpetraExt_BlockMultiVector.h
    block/EpetraExt_BlockVector.h
    block/EpetraExt_BlockUtility.h
    block/EpetraExt_MultiComm.h
    block/EpetraExt_MultiPointModelEvaluator.h
    block/EpetraExt_MultiSerialComm.h
    block/EpetraExt_BlockDiagMatrix.h
    block/EpetraExt_PointToBlockDiagPermute.h
    )

  SET(SOURCES ${SOURCES}
    block/EpetraExt_BlockCrsMatrix.cpp
    block/EpetraExt_BlockMultiVector.cpp
    block/EpetraExt_BlockVector.cpp
    block/EpetraExt_BlockUtility.cpp
    block/EpetraExt_MultiPointModelEvaluator.cpp
    block/EpetraExt_MultiSerialComm.cpp
    block/EpetraExt_BlockDiagMatrix.cpp
    block/EpetraExt_PointToBlockDiagPermute.cpp
    )

  IF(TPL_ENABLE_MPI)

    SET(HEADERS ${HEADERS}
      block/EpetraExt_MultiMpiComm.h
      )

    SET(SOURCES ${SOURCES}
      block/EpetraExt_MultiMpiComm.cpp
      )

  ENDIF()

ENDIF()

#
# src/restrict
#

IF (${PACKAGE_NAME}_BUILD_RESTRICT)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/restrict)

  SET(HEADERS ${HEADERS}
    restrict/EpetraExt_RestrictedMultiVectorWrapper.h
    restrict/EpetraExt_RestrictedCrsMatrixWrapper.h
    )

  SET(SOURCES ${SOURCES}
    restrict/EpetraExt_RestrictedMultiVectorWrapper.cpp
    restrict/EpetraExt_RestrictedCrsMatrixWrapper.cpp
    )

ENDIF()

#
# src/petsc
#

IF (${PACKAGE_NAME}_USING_PETSC)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/petsc)

  SET(HEADERS ${HEADERS}
    petsc/EpetraExt_PETScAIJMatrix.h
    )

  SET(SOURCES ${SOURCES}
    petsc/EpetraExt_PETScAIJMatrix.cpp
    )

ENDIF()

#
# src/hypre
#

IF (TPL_ENABLE_HYPRE)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/hypre)

  SET(HEADERS ${HEADERS}
    hypre/EpetraExt_HypreIJMatrix.h
    )

  SET(SOURCES ${SOURCES}
    hypre/EpetraExt_HypreIJMatrix.cpp
    )

ENDIF()

#
# PyTrilinos support
#

IF (${PROJECT_NAME}_ENABLE_PyTrilinos)

  FIND_PACKAGE(PythonInclude REQUIRED)

  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  epetraext
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
