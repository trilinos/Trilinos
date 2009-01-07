
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
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../)

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

#
# Core headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  fei_Data.hpp
  fei_defs.h
  fei_DirichletBCRecord.hpp
  fei_DirichletBCManager.hpp
  fei_FiniteElementData.hpp
  fei_LinearSystemCore.hpp
  fei_LinSysCore_flexible.hpp
  fei_Lookup.hpp
  fei_BlockDescriptor.hpp
  fei_CSRMat.hpp
  fei_CSVec.hpp
  cfei.h
  fei_ConnectivityTable.hpp
  fei_EqnBuffer.hpp
  fei_EqnCommMgr.hpp
  fei_ErrMacros.hpp
  fei_FEDataFilter.hpp
  FEI.hpp
  feiArray.hpp
  fei_bool.h
  fei_chk_mpi.hpp
  fei_ctg_set.hpp
  fei_EqnComm.hpp
  fei_Factory.hpp
  fei_FillableMat.hpp
  fei_FillableVec.hpp
  FEI_Implementation.hpp
  fei_IndexTable.hpp
  fei_iosfwd.hpp
  fei_iostream.hpp
  fei_LinearProblemManager.hpp
  fei_LinearSystem.hpp
  fei_LogFile.hpp
  fei_Logger.hpp
  fei_LogManager.hpp
  fei_macros.hpp
  fei_Matrix.hpp
  fei_Matrix_core.hpp
  fei_Matrix_Impl.hpp
  fei_Matrix_Local.hpp
  fei_Pool.hpp
  fei_Pool_alloc.hpp
  fei_Vector_Local.hpp
  fei_MatrixReducer.hpp
  fei_MatrixGraph.hpp
  fei_MatrixGraph_Impl2.hpp
  fei_mpi.h
  fei_Param.hpp
  fei_ParameterSet.hpp
  fei_ostream_ops.hpp
  fei_Reducer.hpp
  fei_SharedPtr.hpp
  fei_Solver.hpp
  fei_SparseRowGraph.hpp
  fei_fstream.hpp
  fei_sstream.hpp
  fei_utils.hpp
  fei_Vector.hpp
  fei_Vector_core.hpp
  fei_Vector_Impl.hpp
  fei_VectorReducer.hpp
  fei_VectorSpace.hpp
  fei_version.h
  fei_Filter.hpp
  fei_LibraryWrapper.hpp
  fei_LinSysCoreFilter.hpp
  fei_LinSysCore_struct.h
  fei_NodeCommMgr.hpp
  fei_NodeDatabase.hpp
  fei_NodeDescriptor.hpp
  fei_PatternDescriptor.hpp
  fei_ProcEqns.hpp
  fei_SlaveVariable.hpp
  snl_fei_ArrayUtils.hpp
  snl_fei_BlkSizeMsgHandler.hpp
  snl_fei_BlockMatrixTraits.hpp
  snl_fei_Broker_FEData.hpp
  snl_fei_Broker.hpp
  snl_fei_Broker_LinSysCore.hpp
  fei_CommUtils.hpp
  fei_ConnectivityBlock.hpp
  snl_fei_Constraint.hpp
  snl_fei_Factory.hpp
  snl_fei_FEMatrixTraits_FED.hpp
  snl_fei_FEMatrixTraits.hpp
  snl_fei_FEVectorTraits_FED.hpp
  snl_fei_FEVectorTraits.hpp
  fei_FieldMask.hpp
  fei_Graph.hpp
  fei_GraphReducer.hpp
  fei_Graph_Impl.hpp
  snl_fei_LinearSystem_FEData.hpp
  snl_fei_LinearSystem_General.hpp
  fei_Lookup_Impl.hpp
  snl_fei_MapContig.hpp
  snl_fei_MapTraits.hpp
  snl_fei_MapTraits_specialize.hpp
  snl_fei_MatrixTraits_FEData.hpp
  snl_fei_MatrixTraits.hpp
  snl_fei_MatrixTraits_LinProbMgr.hpp
  snl_fei_MatrixTraits_LinSysCore.hpp
  snl_fei_MatrixTraits_SSMat.hpp
  fei_mpiTraits.hpp
  fei_Pattern.hpp
  snl_fei_PointBlockMap.hpp
  snl_fei_RaggedTable.hpp
  snl_fei_RaggedTable_specialize.hpp
  snl_fei_RecordCollection.hpp
  fei_Record.hpp
  snl_fei_RecordMsgHandler.hpp
  snl_fei_SetTraits.hpp
  snl_fei_SetTraits_specialize.hpp
  fei_SharedIDs.hpp
  SNL_FEI_Structure.hpp
  snl_fei_SubdMsgHandler.hpp
  fei_FEI_Impl.hpp
  fei_TemplateUtils.hpp
  snl_fei_Utils.hpp
  snl_fei_VectorTraits_FEData.hpp
  snl_fei_VectorTraits.hpp
  snl_fei_VectorTraits_LinProbMgr.hpp
  snl_fei_VectorTraits_LinSysCore.hpp
  snl_fei_VectorTraits_SSVec.hpp
  fei_SSGraph.hpp
  fei_SSMat.hpp
  fei_SSVec.hpp
  )

APPEND_SET(SOURCES
  fei_BlockDescriptor.cpp
  fei_CSRMat.cpp
  fei_CSVec.cpp
  fei_DirichletBCManager.cpp
  cfei.cpp
  fei_EqnBuffer.cpp
  fei_EqnCommMgr.cpp
  fei_FEDataFilter.cpp
  FEI_Implementation.cpp
  fei_ostream_ops.cpp
  fei_EqnComm.cpp
  fei_Factory.cpp
  fei_FillableMat.cpp
  fei_FillableVec.cpp
  fei_LinearSystem.cpp
  fei_LogFile.cpp
  fei_Logger.cpp
  fei_LogManager.cpp
  fei_Matrix_core.cpp
  fei_Matrix_Local.cpp
  fei_Vector_Local.cpp
  fei_MatrixReducer.cpp
  fei_MatrixGraph_Impl2.cpp
  fei_Param.cpp
  fei_ParameterSet.cpp
  fei_Pool.cpp
  fei_Pool_alloc.cpp
  fei_Reducer.cpp
  fei_utils.cpp
  fei_VectorSpace.cpp
  fei_Vector_core.cpp
  fei_VectorReducer.cpp
  fei_Filter.cpp
  fei_LibraryWrapper.cpp
  fei_LinSysCoreFilter.cpp
  fei_NodeCommMgr.cpp
  fei_NodeDatabase.cpp
  fei_NodeDescriptor.cpp
  fei_PatternDescriptor.cpp
  fei_ProcEqns.cpp
  fei_Solver.cpp
  snl_fei_BlkSizeMsgHandler.cpp
  snl_fei_Broker_FEData.cpp
  snl_fei_Broker_LinSysCore.cpp
  fei_CommUtils.cpp
  fei_ConnectivityBlock.cpp
  snl_fei_Factory.cpp
  fei_FieldMask.cpp
  fei_GraphReducer.cpp
  fei_Graph_Impl.cpp
  snl_fei_LinearSystem_FEData.cpp
  snl_fei_LinearSystem_General.cpp
  fei_Lookup_Impl.cpp
  fei_Pattern.cpp
  snl_fei_PointBlockMap.cpp
  snl_fei_RaggedTable_specialize.cpp
  fei_Record.cpp
  snl_fei_RecordCollection.cpp
  snl_fei_RecordMsgHandler.cpp
  SNL_FEI_Structure.cpp
  snl_fei_SubdMsgHandler.cpp
  fei_FEI_Impl.cpp
  snl_fei_Utils.cpp
  fei_SSGraph.cpp
  fei_SSMat.cpp
  fei_SSVec.cpp
  )

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  fei
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
