# @HEADER
# ************************************************************************
# 
#              Finite Element Interface to Linear Solvers
#                 Copyright (2001) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# Questions? Contact Alan Williams (william@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside FEI. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package fei subdirectory base built successfully."
	@echo ""

include $(top_builddir)/Makefile.export.fei

DEFS = 
DEFAULT_INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/base -I$(top_builddir) $(TEUCHOS_INCLUDES)

CORE_H = \
	$(top_srcdir)/base/fei_Data.hpp \
	$(top_srcdir)/base/fei_defs.h \
	$(top_srcdir)/base/fei_DirichletBCRecord.hpp \
	$(top_srcdir)/base/fei_DirichletBCManager.hpp \
	$(top_srcdir)/base/fei_FiniteElementData.hpp \
	$(top_srcdir)/base/fei_LinearSystemCore.hpp \
	$(top_srcdir)/base/fei_LinSysCore_flexible.hpp \
	$(top_srcdir)/base/fei_Lookup.hpp \
	$(top_srcdir)/base/fei_BlockDescriptor.hpp \
	$(top_srcdir)/base/fei_CSRMat.hpp \
	$(top_srcdir)/base/fei_CSVec.hpp \
	$(top_srcdir)/base/cfei.h \
	$(top_srcdir)/base/fei_ConnectivityTable.hpp \
	$(top_srcdir)/base/fei_EqnBuffer.hpp \
	$(top_srcdir)/base/fei_EqnCommMgr.hpp \
	$(top_srcdir)/base/fei_ErrMacros.hpp \
	$(top_srcdir)/base/fei_FEDataFilter.hpp \
	$(top_srcdir)/base/FEI.hpp \
	$(top_srcdir)/base/feiArray.hpp \
	$(top_srcdir)/base/fei_bool.h \
	$(top_srcdir)/base/fei_chk_mpi.hpp \
	$(top_srcdir)/base/fei_ctg_set.hpp \
	$(top_srcdir)/base/fei_EqnComm.hpp \
	$(top_srcdir)/base/fei_Factory.hpp \
	$(top_srcdir)/base/fei_FillableMat.hpp \
	$(top_srcdir)/base/fei_FillableVec.hpp \
	$(top_srcdir)/base/FEI_Implementation.hpp \
	$(top_srcdir)/base/fei_IndexTable.hpp \
	$(top_srcdir)/base/fei_iosfwd.hpp \
	$(top_srcdir)/base/fei_iostream.hpp \
	$(top_srcdir)/base/fei_LinearProblemManager.hpp \
	$(top_srcdir)/base/fei_LinearSystem.hpp \
	$(top_srcdir)/base/fei_LogFile.hpp \
	$(top_srcdir)/base/fei_Logger.hpp \
	$(top_srcdir)/base/fei_LogManager.hpp \
	$(top_srcdir)/base/fei_macros.hpp \
	$(top_srcdir)/base/fei_Matrix.hpp \
	$(top_srcdir)/base/fei_Matrix_core.hpp \
	$(top_srcdir)/base/fei_Matrix_Impl.hpp \
	$(top_srcdir)/base/fei_Matrix_Local.hpp \
	$(top_srcdir)/base/fei_Pool.hpp \
	$(top_srcdir)/base/fei_Pool_alloc.hpp \
	$(top_srcdir)/base/fei_Vector_Local.hpp \
	$(top_srcdir)/base/fei_MatrixReducer.hpp \
	$(top_srcdir)/base/fei_MatrixGraph.hpp \
	$(top_srcdir)/base/fei_MatrixGraph_Impl2.hpp \
	$(top_srcdir)/base/fei_mpi.h \
	$(top_srcdir)/base/fei_Param.hpp \
	$(top_srcdir)/base/fei_ParameterSet.hpp \
	$(top_srcdir)/base/fei_ostream_ops.hpp \
	$(top_srcdir)/base/fei_Reducer.hpp \
	$(top_srcdir)/base/fei_SharedPtr.hpp \
	$(top_srcdir)/base/fei_Solver.hpp \
	$(top_srcdir)/base/fei_SparseRowGraph.hpp \
	$(top_srcdir)/base/fei_fstream.hpp \
	$(top_srcdir)/base/fei_sstream.hpp \
	$(top_srcdir)/base/fei_utils.hpp \
	$(top_srcdir)/base/fei_Vector.hpp \
	$(top_srcdir)/base/fei_Vector_core.hpp \
	$(top_srcdir)/base/fei_Vector_Impl.hpp \
	$(top_srcdir)/base/fei_VectorReducer.hpp \
	$(top_srcdir)/base/fei_VectorSpace.hpp \
	$(top_srcdir)/base/fei_version.h \
	$(top_srcdir)/base/fei_Filter.hpp \
	$(top_srcdir)/base/fei_LibraryWrapper.hpp \
	$(top_srcdir)/base/fei_LinSysCoreFilter.hpp \
	$(top_srcdir)/base/fei_LinSysCore_struct.h \
	$(top_srcdir)/base/fei_NodeCommMgr.hpp \
	$(top_srcdir)/base/fei_NodeDatabase.hpp \
	$(top_srcdir)/base/fei_NodeDescriptor.hpp \
	$(top_srcdir)/base/fei_PatternDescriptor.hpp \
	$(top_srcdir)/base/fei_ProcEqns.hpp \
	$(top_srcdir)/base/fei_SlaveVariable.hpp \
	$(top_srcdir)/base/snl_fei_ArrayUtils.hpp \
	$(top_srcdir)/base/snl_fei_BlkSizeMsgHandler.hpp \
	$(top_srcdir)/base/snl_fei_BlockMatrixTraits.hpp \
	$(top_srcdir)/base/snl_fei_Broker_FEData.hpp \
	$(top_srcdir)/base/snl_fei_Broker.hpp \
	$(top_srcdir)/base/snl_fei_Broker_LinSysCore.hpp \
	$(top_srcdir)/base/fei_CommUtils.hpp \
	$(top_srcdir)/base/fei_ConnectivityBlock.hpp \
	$(top_srcdir)/base/snl_fei_Constraint.hpp \
	$(top_srcdir)/base/snl_fei_Factory.hpp \
	$(top_srcdir)/base/snl_fei_FEMatrixTraits_FED.hpp \
	$(top_srcdir)/base/snl_fei_FEMatrixTraits.hpp \
	$(top_srcdir)/base/snl_fei_FEVectorTraits_FED.hpp \
	$(top_srcdir)/base/snl_fei_FEVectorTraits.hpp \
	$(top_srcdir)/base/fei_FieldMask.hpp \
	$(top_srcdir)/base/fei_Graph.hpp \
	$(top_srcdir)/base/fei_GraphReducer.hpp \
	$(top_srcdir)/base/fei_Graph_Impl.hpp \
	$(top_srcdir)/base/snl_fei_LinearSystem_FEData.hpp \
	$(top_srcdir)/base/snl_fei_LinearSystem_General.hpp \
	$(top_srcdir)/base/fei_Lookup_Impl.hpp \
	$(top_srcdir)/base/snl_fei_MapContig.hpp \
	$(top_srcdir)/base/snl_fei_MapTraits.hpp \
	$(top_srcdir)/base/snl_fei_MapTraits_specialize.hpp \
	$(top_srcdir)/base/snl_fei_MatrixTraits_FEData.hpp \
	$(top_srcdir)/base/snl_fei_MatrixTraits.hpp \
	$(top_srcdir)/base/snl_fei_MatrixTraits_LinProbMgr.hpp \
	$(top_srcdir)/base/snl_fei_MatrixTraits_LinSysCore.hpp \
	$(top_srcdir)/base/snl_fei_MatrixTraits_SSMat.hpp \
	$(top_srcdir)/base/fei_mpiTraits.hpp \
	$(top_srcdir)/base/fei_Pattern.hpp \
	$(top_srcdir)/base/snl_fei_PointBlockMap.hpp \
	$(top_srcdir)/base/snl_fei_RaggedTable.hpp \
	$(top_srcdir)/base/snl_fei_RaggedTable_specialize.hpp \
	$(top_srcdir)/base/snl_fei_RecordCollection.hpp \
	$(top_srcdir)/base/fei_Record.hpp \
	$(top_srcdir)/base/snl_fei_RecordMsgHandler.hpp \
	$(top_srcdir)/base/snl_fei_SetTraits.hpp \
	$(top_srcdir)/base/snl_fei_SetTraits_specialize.hpp \
	$(top_srcdir)/base/fei_SharedIDs.hpp \
	$(top_srcdir)/base/SNL_FEI_Structure.hpp \
	$(top_srcdir)/base/snl_fei_SubdMsgHandler.hpp \
	$(top_srcdir)/base/fei_FEI_Impl.hpp \
	$(top_srcdir)/base/fei_TemplateUtils.hpp \
	$(top_srcdir)/base/snl_fei_Utils.hpp \
	$(top_srcdir)/base/snl_fei_VectorTraits_FEData.hpp \
	$(top_srcdir)/base/snl_fei_VectorTraits.hpp \
	$(top_srcdir)/base/snl_fei_VectorTraits_LinProbMgr.hpp \
	$(top_srcdir)/base/snl_fei_VectorTraits_LinSysCore.hpp \
	$(top_srcdir)/base/snl_fei_VectorTraits_SSVec.hpp \
	$(top_srcdir)/base/fei_SSGraph.hpp \
	$(top_srcdir)/base/fei_SSMat.hpp \
	$(top_srcdir)/base/fei_SSVec.hpp

CORE = \
	$(srcdir)/fei_BlockDescriptor.cpp \
	$(srcdir)/fei_CSRMat.cpp \
	$(srcdir)/fei_CSVec.cpp \
	$(srcdir)/fei_DirichletBCManager.cpp \
	$(srcdir)/cfei.cpp \
	$(srcdir)/fei_EqnBuffer.cpp \
	$(srcdir)/fei_EqnCommMgr.cpp \
	$(srcdir)/fei_FEDataFilter.cpp \
	$(srcdir)/FEI_Implementation.cpp \
	$(srcdir)/fei_ostream_ops.cpp \
	$(srcdir)/fei_EqnComm.cpp \
	$(srcdir)/fei_Factory.cpp \
	$(srcdir)/fei_FillableMat.cpp \
	$(srcdir)/fei_FillableVec.cpp \
	$(srcdir)/fei_LinearSystem.cpp \
	$(srcdir)/fei_LogFile.cpp \
	$(srcdir)/fei_Logger.cpp \
	$(srcdir)/fei_LogManager.cpp \
	$(srcdir)/fei_Matrix_core.cpp \
	$(srcdir)/fei_Matrix_Local.cpp \
	$(srcdir)/fei_Vector_Local.cpp \
	$(srcdir)/fei_MatrixReducer.cpp \
	$(srcdir)/fei_MatrixGraph_Impl2.cpp \
	$(srcdir)/fei_Param.cpp \
	$(srcdir)/fei_ParameterSet.cpp \
	$(srcdir)/fei_Pool.cpp \
	$(srcdir)/fei_Pool_alloc.cpp \
	$(srcdir)/fei_Reducer.cpp \
	$(srcdir)/fei_utils.cpp \
	$(srcdir)/fei_VectorSpace.cpp \
	$(srcdir)/fei_Vector_core.cpp \
	$(srcdir)/fei_VectorReducer.cpp \
	$(srcdir)/fei_Filter.cpp \
	$(srcdir)/fei_LibraryWrapper.cpp \
	$(srcdir)/fei_LinSysCoreFilter.cpp \
	$(srcdir)/fei_NodeCommMgr.cpp \
	$(srcdir)/fei_NodeDatabase.cpp \
	$(srcdir)/fei_NodeDescriptor.cpp \
	$(srcdir)/fei_PatternDescriptor.cpp \
	$(srcdir)/fei_ProcEqns.cpp \
	$(srcdir)/fei_Solver.cpp \
	$(srcdir)/snl_fei_BlkSizeMsgHandler.cpp \
	$(srcdir)/snl_fei_Broker_FEData.cpp \
	$(srcdir)/snl_fei_Broker_LinSysCore.cpp \
	$(srcdir)/fei_CommUtils.cpp \
	$(srcdir)/fei_ConnectivityBlock.cpp \
	$(srcdir)/snl_fei_Factory.cpp \
	$(srcdir)/fei_FieldMask.cpp \
	$(srcdir)/fei_GraphReducer.cpp \
	$(srcdir)/fei_Graph_Impl.cpp \
	$(srcdir)/snl_fei_LinearSystem_FEData.cpp \
	$(srcdir)/snl_fei_LinearSystem_General.cpp \
	$(srcdir)/fei_Lookup_Impl.cpp \
	$(srcdir)/fei_Pattern.cpp \
	$(srcdir)/snl_fei_PointBlockMap.cpp \
	$(srcdir)/snl_fei_RaggedTable_specialize.cpp \
	$(srcdir)/fei_Record.cpp \
	$(srcdir)/snl_fei_RecordCollection.cpp \
	$(srcdir)/snl_fei_RecordMsgHandler.cpp \
	$(srcdir)/SNL_FEI_Structure.cpp \
	$(srcdir)/snl_fei_SubdMsgHandler.cpp \
	$(srcdir)/fei_FEI_Impl.cpp \
	$(srcdir)/snl_fei_Utils.cpp \
	$(srcdir)/fei_SSGraph.cpp \
	$(srcdir)/fei_SSMat.cpp \
	$(srcdir)/fei_SSVec.cpp

# ------------------------------------------------------------------------
# FEI library specifications
# ------------------------------------------------------------------------

FEI_LIB = libfei_base.a


FEI_H = \
	$(CORE_H)

libfei_base_a_SOURCES = \
	$(CORE) 

EXTRA_libfei_base_a_SOURCES = 

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR if USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libfei_base_a_AR = $(ALTERNATE_AR)
else

libfei_base_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# On alpha-dec machines, the extra .o's are needed for Epetra only.
# ------------------------------------------------------------------------

libfei_base_a_LIBADD = $(FEI_TEMPLATE_OBJS) $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(FEI_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(FEI_H) 

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

