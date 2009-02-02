
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
../fei_fwd.hpp
../fei_base.hpp
cfei.h
feiArray.hpp
fei_BlockDescriptor.hpp
fei_BlockLinearProblemManager.hpp
fei_bool.h
fei_chk_mpi.hpp
fei_CommUtils.hpp
fei_ConnectivityBlock.hpp
fei_ConnectivityTable.hpp
fei_CSRMat.hpp
fei_CSVec.hpp
fei_ctg_set.hpp
fei_Data.hpp
fei_defs.h
fei_DirichletBCManager.hpp
fei_DirichletBCRecord.hpp
fei_EqnBuffer.hpp
fei_EqnComm.hpp
fei_EqnCommMgr.hpp
fei_EqnRecord.hpp
fei_ErrMacros.hpp
fei_Factory.hpp
fei_FEDataFilter.hpp
fei_FEI_Impl.hpp
fei_FieldMask.hpp
fei_FillableMat.hpp
fei_FillableVec.hpp
fei_Filter.hpp
fei_FiniteElementData.hpp
fei_fstream.hpp
fei_Graph.hpp
fei_Graph_Impl.hpp
fei_GraphReducer.hpp
FEI.hpp
FEI_Implementation.hpp
fei_impl_utils.hpp
fei_IndexTable.hpp
fei_iosfwd.hpp
fei_iostream.hpp
fei_LibraryWrapper.hpp
fei_LinearProblemManager.hpp
fei_LinearSystemCore.hpp
fei_LinearSystem.hpp
fei_LinSysCoreFilter.hpp
fei_LinSysCore_flexible.hpp
fei_LinSysCore_struct.h
fei_LogFile.hpp
fei_Logger.hpp
fei_LogManager.hpp
fei_Lookup.hpp
fei_Lookup_Impl.hpp
fei_macros.hpp
fei_Matrix_core.hpp
fei_MatrixGraph.hpp
fei_MatrixGraph_Impl2.hpp
fei_Matrix.hpp
fei_Matrix_Impl.hpp
fei_Matrix_Local.hpp
fei_MatrixReducer.hpp
fei_MatrixTraits_FEData.hpp
fei_MatrixTraits_FillableMat.hpp
fei_MatrixTraits.hpp
fei_MatrixTraits_LinProbMgr.hpp
fei_MatrixTraits_LinSysCore.hpp
fei_mpi.h
fei_mpiTraits.hpp
fei_NodeCommMgr.hpp
fei_NodeDatabase.hpp
fei_NodeDescriptor.hpp
fei_ostream_ops.hpp
fei_ParameterSet.hpp
fei_Param.hpp
fei_PatternDescriptor.hpp
fei_Pattern.hpp
fei_Pool_alloc.hpp
fei_Pool.hpp
fei_ProcEqns.hpp
fei_Record.hpp
fei_Reducer.hpp
fei_ReverseMapper.hpp
fei_SharedIDs.hpp
fei_SharedPtr.hpp
fei_SlaveVariable.hpp
fei_Solver.hpp
fei_SparseRowGraph.hpp
fei_SSGraph.hpp
fei_sstream.hpp
fei_SSVec.hpp
fei_TemplateUtils.hpp
fei_utils.hpp
fei_Vector_core.hpp
fei_Vector.hpp
fei_Vector_Impl.hpp
fei_Vector_Local.hpp
fei_VectorReducer.hpp
fei_VectorSpace.hpp
fei_VectorTraits_FEData.hpp
fei_VectorTraits_FillableVec.hpp
fei_VectorTraits.hpp
fei_VectorTraits_LinProbMgr.hpp
fei_VectorTraits_LinSysCore.hpp
fei_VectorTraits_SSVec.hpp
fei_version.h
snl_fei_ArrayUtils.hpp
snl_fei_BlkSizeMsgHandler.hpp
snl_fei_BlockMatrixTraits.hpp
snl_fei_Broker_FEData.hpp
snl_fei_Broker.hpp
snl_fei_Broker_LinSysCore.hpp
snl_fei_Constraint.hpp
snl_fei_Factory.hpp
snl_fei_FEMatrixTraits_FED.hpp
snl_fei_FEMatrixTraits.hpp
snl_fei_FEVectorTraits_FED.hpp
snl_fei_FEVectorTraits.hpp
snl_fei_LinearSystem_FEData.hpp
snl_fei_LinearSystem_General.hpp
snl_fei_MapContig.hpp
snl_fei_MapTraits.hpp
snl_fei_MapTraits_specialize.hpp
snl_fei_PointBlockMap.hpp
snl_fei_RaggedTable.hpp
snl_fei_RaggedTable_specialize.hpp
snl_fei_RecordCollection.hpp
snl_fei_RecordMsgHandler.hpp
snl_fei_SetTraits.hpp
snl_fei_SetTraits_specialize.hpp
SNL_FEI_Structure.hpp
snl_fei_SubdMsgHandler.hpp
snl_fei_Utils.hpp
  )

APPEND_SET(SOURCES
cfei.cpp
fei_BlockDescriptor.cpp
fei_CommUtils.cpp
fei_ConnectivityBlock.cpp
fei_CSRMat.cpp
fei_CSVec.cpp
fei_DirichletBCManager.cpp
fei_EqnBuffer.cpp
fei_EqnComm.cpp
fei_EqnCommMgr.cpp
fei_Factory.cpp
fei_FEDataFilter.cpp
fei_FEI_Impl.cpp
fei_FieldMask.cpp
fei_FillableMat.cpp
fei_FillableVec.cpp
fei_Filter.cpp
fei_Graph_Impl.cpp
fei_GraphReducer.cpp
FEI_Implementation.cpp
fei_impl_utils.cpp
fei_LibraryWrapper.cpp
fei_LinearSystem.cpp
fei_LinSysCoreFilter.cpp
fei_LogFile.cpp
fei_Logger.cpp
fei_LogManager.cpp
fei_Lookup_Impl.cpp
fei_Matrix_core.cpp
fei_MatrixGraph_Impl2.cpp
fei_Matrix_Local.cpp
fei_MatrixReducer.cpp
fei_NodeCommMgr.cpp
fei_NodeDatabase.cpp
fei_NodeDescriptor.cpp
fei_ostream_ops.cpp
fei_Param.cpp
fei_ParameterSet.cpp
fei_Pattern.cpp
fei_PatternDescriptor.cpp
fei_Pool_alloc.cpp
fei_Pool.cpp
fei_ProcEqns.cpp
fei_Record.cpp
fei_Reducer.cpp
fei_Solver.cpp
fei_SSGraph.cpp
fei_SSVec.cpp
fei_utils.cpp
fei_Vector_core.cpp
fei_Vector_Local.cpp
fei_VectorReducer.cpp
fei_VectorSpace.cpp
snl_fei_BlkSizeMsgHandler.cpp
snl_fei_Broker_FEData.cpp
snl_fei_Broker_LinSysCore.cpp
snl_fei_Factory.cpp
snl_fei_LinearSystem_FEData.cpp
snl_fei_LinearSystem_General.cpp
snl_fei_PointBlockMap.cpp
snl_fei_RaggedTable_specialize.cpp
snl_fei_RecordCollection.cpp
snl_fei_RecordMsgHandler.cpp
SNL_FEI_Structure.cpp
snl_fei_SubdMsgHandler.cpp
snl_fei_Utils.cpp
  )

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  fei_base
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
