# @HEADER
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package ml subdirectory src built successfully."
	@echo ""

#np# Note the following paragraph.
# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside ml. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
#np# There is, at present, no distinction between most of the categories, 
#np# other than those which start with EXTRA_.  All categories which do not
#np# start with EXTRA_ are included in all builds.  Categories which start
#np# with EXTRA_ are included in some builds (typically controlled by 
#np# --with-xyz or --without-xyz)
#
#np# Make sure to list all source files in one of the following categories.
CORE_H = \
$(srcdir)/Include/ml_defs.h     \
$(srcdir)/Comm/ml_comm.h \
$(srcdir)/Include/ml_include.h \
$(srcdir)/Include/ml_common.h \
$(srcdir)/Coarsen/ml_amg_genP.h \
$(srcdir)/Coarsen/ml_amg.h \
$(srcdir)/Coarsen/ml_agg_genP.h \
$(srcdir)/Coarsen/ml_agg_min_energy.h \
$(srcdir)/Coarsen/ml_check.h \
$(srcdir)/Coarsen/ml_aggregate.h \
$(srcdir)/Coarsen/ml_ggraph.h \
$(srcdir)/FEGrid/ml_mapper.h \
$(srcdir)/FEGrid/ml_gridagx.h \
$(srcdir)/FEGrid/ml_pde.h       $(srcdir)/FEGrid/ml_bdrypts.h \
$(srcdir)/Coarsen/ml_agg_reitzinger.h \
$(srcdir)/Coarsen/ml_agg_METIS.h \
$(srcdir)/Coarsen/ml_agg_user.h \
$(srcdir)/Coarsen/ml_agg_VBMETIS.h \
$(srcdir)/Coarsen/ml_agg_Zoltan.h \
$(srcdir)/Coarsen/ml_agg_ParMETIS.h \
$(srcdir)/Coarsen/ml_qr_fix.h       \
$(srcdir)/FEGrid/ml_gridfunc.h \
$(srcdir)/FEGrid/ml_setup.h \
$(srcdir)/FEGrid/ml_elementagx.h  $(srcdir)/FEGrid/ml_grid.h \
$(srcdir)/Main/ml_struct.h $(srcdir)/Main/ml_seg_precond.h \
$(srcdir)/Main/ml_1level.h \
$(srcdir)/Main/mli_solver.h \
$(srcdir)/Main/ml_ggb.h \
$(srcdir)/Main/ml_viz_stats.h \
$(srcdir)/Comm/ml_comminfoagx.h \
$(srcdir)/Comm/ml_comminfoop.h \
$(srcdir)/Smoother/ml_csolve.h \
$(srcdir)/Smoother/ml_mls.h \
$(srcdir)/Smoother/ml_smoother.h \
$(srcdir)/Smoother/ml_solver.h \
$(srcdir)/Smoother/ml_superlu.h \
$(srcdir)/Smoother/ml_xxt.h \
$(srcdir)/Smoother/ml_xyt.h \
$(srcdir)/Krylov/ml_bicgstabl.h \
$(srcdir)/Krylov/ml_cg.h \
$(srcdir)/Krylov/ml_gmres.h \
$(srcdir)/Krylov/ml_krylov.h \
$(srcdir)/Operator/ml_mat_formats.h \
$(srcdir)/Operator/ml_operatoragx.h \
$(srcdir)/Operator/ml_operator_blockmat.h \
$(srcdir)/Operator/ml_operator.h \
$(srcdir)/Operator/ml_op_utils.h \
$(srcdir)/Operator/ml_rap.h

CORE = \
$(srcdir)/Coarsen/ml_agg_uncoupled.c \
$(srcdir)/Coarsen/ml_amg.c \
$(srcdir)/Coarsen/ml_agg_coupled.c     $(srcdir)/Coarsen/ml_amg_genP.c \
$(srcdir)/Coarsen/ml_check.c \
$(srcdir)/Coarsen/ml_agg_dd.c \
$(srcdir)/Coarsen/ml_agg_genP.c \
$(srcdir)/Coarsen/ml_agg_min_energy.cpp \
$(srcdir)/Coarsen/ml_amg_MIS.c \
$(srcdir)/Coarsen/ml_agg_MIS.c \
$(srcdir)/Coarsen/ml_agg_METIS.c \
$(srcdir)/Coarsen/ml_agg_user.c \
$(srcdir)/Coarsen/ml_agg_VBMETIS.c \
$(srcdir)/Coarsen/ml_agg_Zoltan.c \
$(srcdir)/Coarsen/ml_agg_ParMETIS.c \
$(srcdir)/Coarsen/ml_aggregate.c \
$(srcdir)/Coarsen/ml_ggraph.c \
$(srcdir)/Coarsen/ml_agg_reitzinger.c  \
$(srcdir)/Coarsen/ml_qr_fix.c \
 \
$(srcdir)/FEGrid/ml_gridagx.c \
$(srcdir)/FEGrid/ml_pde.c \
$(srcdir)/FEGrid/ml_bdrypts.c   $(srcdir)/FEGrid/ml_grid.c \
$(srcdir)/FEGrid/ml_gridfunc.c  $(srcdir)/FEGrid/ml_setup.c \
$(srcdir)/FEGrid/ml_elementagx.c \
$(srcdir)/FEGrid/ml_get_basis.c   $(srcdir)/FEGrid/ml_mapper.c \
 \
$(srcdir)/Main/ml_struct.c $(srcdir)/Main/driver.c \
$(srcdir)/Main/mli_solver.c \
$(srcdir)/Main/ml_seg_precond.c \
$(srcdir)/Main/ml_ggb.c \
 \
$(srcdir)/Comm/ml_comm.c \
$(srcdir)/Comm/ml_comminfoagx.c \
$(srcdir)/Comm/ml_comminfoop.c \
$(srcdir)/Comm/ml_exch_row.c \
 \
$(srcdir)/Smoother/ml_csolve.c \
$(srcdir)/Smoother/ml_smoother.c \
$(srcdir)/Smoother/ml_solver.c \
$(srcdir)/Smoother/ml_superlu.c \
$(srcdir)/Smoother/ml_xxt.c \
$(srcdir)/Smoother/ml_xyt.c \
 \
$(srcdir)/Krylov/ml_bicgstabl.c \
$(srcdir)/Krylov/ml_cg.c \
$(srcdir)/Krylov/ml_gmres.c \
$(srcdir)/Krylov/ml_krylov.c \
 \
$(srcdir)/Operator/ml_matmat_mult.c \
$(srcdir)/Operator/ml_operatoragx.c \
$(srcdir)/Operator/ml_operator_blockmat.c \
$(srcdir)/Operator/ml_operator.c \
$(srcdir)/Operator/ml_op_utils.c \
$(srcdir)/Operator/ml_rap.c \
$(srcdir)/Operator/ml_mat_formats.c $(srcdir)/Operator/ml_rap_utils.c

UTIL_H = \
$(srcdir)/Utils/ml_aztec_utils.h      \
$(srcdir)/Utils/ml_epetra_utils.h     \
$(srcdir)/Utils/ml_utils.h            \
$(srcdir)/Utils/ml_memory.h           \
$(srcdir)/Utils/ml_agg_info.h         \
$(srcdir)/Utils/ml_aztec_lapack.h     \
$(srcdir)/Utils/ml_rbm.h              \
$(srcdir)/Utils/ml_read_utils.h       \
$(srcdir)/Utils/ml_superlu_lapack.h   \
$(srcdir)/Utils/ml_epetra_operator.h  \
$(srcdir)/Utils/ml_twogrid_analysis.h \
$(srcdir)/Utils/ml_intlist.h          \
$(srcdir)/Utils/ml_vec.h    \
$(srcdir)/Utils/ml_lapack.h \
$(srcdir)/Utils/ml_utils.h  \
$(srcdir)/Utils/ml_eigf2c.h \
$(srcdir)/Utils/ml_vampir.h \
$(srcdir)/Utils/ml_viz_opendx.h   \
$(srcdir)/Utils/ml_viz_xyz.h      \
$(srcdir)/Utils/ml_viz_vtk.h      \
$(srcdir)/Utils/ml_amesos.h       \
$(srcdir)/Utils/ml_amesos_wrap.h  \
$(srcdir)/Utils/ml_ifpack.h       \
$(srcdir)/Utils/ml_ifpack_wrap.h  \
$(srcdir)/Utils/ml_self.h         \
$(srcdir)/Utils/ml_self_wrap.h    \
$(srcdir)/Utils/ml_epetra_preconditioner.h               \
$(srcdir)/Utils/ml_anasazi.h                             \
$(srcdir)/Utils/ml_epetra.h                              \
$(srcdir)/Utils/ml_MultiLevelOperator.h                  \
$(srcdir)/Utils/ml_MultiLevelPreconditioner.h            \
$(srcdir)/Utils/ml_RowMatrix.h                           \
$(srcdir)/Utils/ml_FilterType.h                          \
$(srcdir)/Utils/ml_DD_prec.h                             \
$(srcdir)/Utils/ml_Ifpack_ML.h                           \
$(srcdir)/Utils/ml_Epetra_wrap_CrsGraph_as_RowMatrix.H   \
$(srcdir)/Utils/ml_Preconditioner.h                      \
$(srcdir)/Utils/Epetra_Operator_With_MatMat.h            \
$(srcdir)/Utils/Epetra_Multi_CrsMatrix.h                 \
$(srcdir)/RefMaxwell/ml_EdgeMatrixFreePreconditioner.h   \
$(srcdir)/RefMaxwell/ml_RefMaxwell_11_Operator.h         \
$(srcdir)/RefMaxwell/ml_RefMaxwell.h                     \
$(srcdir)/MatrixFree/ml_MatrixFreePreconditioner.h       \
$(srcdir)/MatrixFree/ml_ElementByElement_SingleElement.h \
$(srcdir)/NonlinML/nlnml_ConstrainedMultiLevelOperator.H \
$(srcdir)/NonlinML/nlnml_finelevelnoxinterface.H         \
$(srcdir)/NonlinML/nlnml_coarselevelnoxinterface.H       \
$(srcdir)/NonlinML/nlnml_preconditioner.H                \
$(srcdir)/NonlinML/nlnml_preconditioner_utils.H          \
$(srcdir)/NonlinML/nlnml_linearsystem.H                  \
$(srcdir)/NonlinML/nlnml_nonlinearlevel.H                \
$(srcdir)/NonlinML/nlnml_prepostoperator.H               \
$(srcdir)/MLAPI/MLAPI.h                      \
$(srcdir)/MLAPI/MLAPI_Aggregation.h          \
$(srcdir)/MLAPI/MLAPI_BaseObject.h           \
$(srcdir)/MLAPI/MLAPI_BaseOperator.h         \
$(srcdir)/MLAPI/MLAPI_CompObject.h           \
$(srcdir)/MLAPI/MLAPI_Defaults.h             \
$(srcdir)/MLAPI/MLAPI_Eig.h                  \
$(srcdir)/MLAPI/MLAPI_EpetraBaseOperator.h   \
$(srcdir)/MLAPI/MLAPI_Expressions.h          \
$(srcdir)/MLAPI/MLAPI_Gallery.h              \
$(srcdir)/MLAPI/MLAPI_InverseOperator.h      \
$(srcdir)/MLAPI/MLAPI_Krylov.h               \
$(srcdir)/MLAPI/MLAPI_MATLABStream.h         \
$(srcdir)/MLAPI/MLAPI_SerialMatrix.h         \
$(srcdir)/MLAPI/MLAPI_DistributedMatrix.h    \
$(srcdir)/MLAPI/MLAPI_MultiLevelSA.h         \
$(srcdir)/MLAPI/MLAPI_MultiLevelAdaptiveSA.h \
$(srcdir)/MLAPI/MLAPI_MultiVector.h          \
$(srcdir)/MLAPI/MLAPI_MultiVector_Utils.h    \
$(srcdir)/MLAPI/MLAPI_Operator_Box.h         \
$(srcdir)/MLAPI/MLAPI_Operator.h             \
$(srcdir)/MLAPI/MLAPI_Operator_Utils.h       \
$(srcdir)/MLAPI/MLAPI_Space.h                \
$(srcdir)/MLAPI/MLAPI_TimeObject.h           \
$(srcdir)/MLAPI/MLAPI_Workspace.h            \
$(srcdir)/MLAPI/MLAPI_Error.h                \
$(srcdir)/MLAPI/MLAPI_SAMIS.h                \
$(srcdir)/MLAPI/MLAPI_LinearCombinations.h   \
$(srcdir)/MLAPI/MLAPI_BaseLinearCombination.h

UTIL = \
$(srcdir)/Utils/ml_rbm.c                                         \
$(srcdir)/Utils/ml_aztec_utils.c                                 \
$(srcdir)/Utils/ml_read_utils.c                                  \
$(srcdir)/Utils/ml_twogrid_analysis.c                            \
$(srcdir)/Utils/ml_epetra_utils.cpp                              \
$(srcdir)/Utils/ml_utils.c                                       \
$(srcdir)/Utils/ml_intlist.c                                     \
$(srcdir)/Utils/ml_vec.c                                         \
$(srcdir)/Utils/ml_memory.c                                      \
$(srcdir)/Utils/ml_eigf2c.c                                      \
$(srcdir)/Utils/ml_agg_info.c                                    \
$(srcdir)/Utils/ml_viz_opendx.c                                  \
$(srcdir)/Utils/ml_viz_xyz.c                                     \
$(srcdir)/Utils/ml_viz_vtk.c                                     \
$(srcdir)/Utils/ml_amesos.c                                      \
$(srcdir)/Utils/ml_amesos_wrap.cpp                               \
$(srcdir)/Utils/ml_ifpack.c                                      \
$(srcdir)/Utils/ml_ifpack_wrap.cpp                               \
$(srcdir)/Utils/ml_self.c                                        \
$(srcdir)/Utils/ml_self_wrap.cpp                                 \
$(srcdir)/Utils/ml_MultiLevelPreconditioner.cpp                  \
$(srcdir)/Utils/ml_MultiLevelPreconditioner_Filtering.cpp        \
$(srcdir)/Utils/ml_MultiLevelPreconditioner_NullSpace.cpp        \
$(srcdir)/Utils/ml_MultiLevelPreconditioner_Defaults.cpp         \
$(srcdir)/Utils/ml_MultiLevelPreconditioner_Smoothers.cpp        \
$(srcdir)/Utils/ml_MultiLevelPreconditioner_Viz.cpp              \
$(srcdir)/Utils/ml_MultiLevelPreconditioner_Adapt.cpp            \
$(srcdir)/Utils/ml_MultiLevelPreconditioner_Analyze.cpp          \
$(srcdir)/Utils/ml_MultiLevelPreconditioner_Test.cpp             \
$(srcdir)/Utils/ml_MultiLevelPreconditioner_Aux.cpp              \
$(srcdir)/Utils/ml_MultiLevelPreconditioner_XML.cpp              \
$(srcdir)/Utils/ml_anasazi.cpp                                   \
$(srcdir)/Utils/ml_MultiLevelOperator.cpp                        \
$(srcdir)/Utils/ml_RowMatrix.cpp                                 \
$(srcdir)/Utils/ml_DD_prec.cpp                                   \
$(srcdir)/Utils/Epetra_Multi_CrsMatrix.cpp                       \
$(srcdir)/RefMaxwell/ml_EdgeMatrixFreePreconditioner.cpp         \
$(srcdir)/RefMaxwell/ml_RefMaxwell_11_Operator.cpp               \
$(srcdir)/RefMaxwell/ml_RefMaxwell.cpp                           \
$(srcdir)/MatrixFree/ml_MatrixFreePreconditioner.cpp             \
$(srcdir)/MLAPI/MLAPI_Defaults.cpp                               \
$(srcdir)/MLAPI/MLAPI_Gallery.cpp                                \
$(srcdir)/MLAPI/MLAPI_InverseOperator.cpp                        \
$(srcdir)/MLAPI/MLAPI_Operator_Utils.cpp                         \
$(srcdir)/MLAPI/MLAPI_Eig.cpp                                    \
$(srcdir)/MLAPI/MLAPI_Krylov.cpp                                 \
$(srcdir)/MLAPI/MLAPI_Workspace.cpp                              \
$(srcdir)/MLAPI/MLAPI_Expressions.cpp                            \
$(srcdir)/MLAPI/MLAPI_MultiVector_Utils.cpp                      \
$(srcdir)/MLAPI/MLAPI_BaseObject.cpp                             \
$(srcdir)/MLAPI/MLAPI_Aggregation.cpp                            \
$(srcdir)/MLAPI/MLAPI_Error.cpp                                  \
$(srcdir)/MLAPI/MLAPI_DistributedMatrix.cpp                      \
$(srcdir)/MLAPI/MLAPI_SAMIS.cpp                                  \
$(srcdir)/MLAPI/MLAPI_LinearCombinations.cpp                     \
$(srcdir)/NonlinML/nlnml_prepostoperator.cpp                     \
$(srcdir)/NonlinML/nlnml_ConstrainedMultiLevelOperator.cpp       \
$(srcdir)/NonlinML/nlnml_linearsystem.cpp                        \
$(srcdir)/NonlinML/nlnml_nonlinearlevel.cpp                      \
$(srcdir)/NonlinML/nlnml_coarselevelnoxinterface.cpp             \
$(srcdir)/NonlinML/nlnml_preconditioner1.cpp                     \
$(srcdir)/NonlinML/nlnml_preconditioner_utils.cpp

ENUM_H = 

UNIPROC_H = 

UNIPROC = 

SERIAL_H = 

SERIAL = 

#
#np#  Files which are included in select builds are typically have
#np#  names that start with EXTRA
#
#np#  EXTRA_MPI is control by --with-mpi and --without-mpi which is
#np#  set up by TAC_ARC_CONFIG_MPI 
#
EXTRA_MPI_H =

EXTRA_MPI = 

if HAVE_MPI 
MPI_H = \
	$(EXTRA_MPI_H)
MPI = \
	$(EXTRA_MPI)

else

MPI_H =
MPI =

endif
#np# This section should come out, or be modified to reflect options in the 
#np# package being autotool'ed
#EXTRA_NEWP_SWAHILI_H = \
#	$(srcdir)/Newp_Jambo.h 
#
#EXTRA_NEWP_SWAHILI = \
#	$(srcdir)/Newp_Jambo.cpp
#
#if HAVE_NEWP_SWAHILI
#
#NEWP_SWAHILI_H = \
#	$(EXTRA_NEWP_SWAHILI_H)
#NEWP_SWAHILI = \
#	$(EXTRA_NEWP_SWAHILI)
#
#else
#
#NEWP_SWAHILI_H =
#NEWP_SWAHILI =
#
#endif



USERBASE_H = 


LANG_WRAPS_H = 

LANG_WRAPS = 

DIST_H = 

DIST = 

COMP_H =

COMP = 


#np# List locations of include files within Trilinos but outside of this 
#np# package
#np# - Add packages that this package is dependent on
#np#   Remove epetra if this package does not need it

if ML_USING_NOX
nox_incdir = -I$(top_srcdir)/../nox/src -I$(top_srcdir)/../nox/src-epetra -I$(top_builddir)/../nox/src
else
nox_incdir =
endif

include $(top_builddir)/Makefile.export.ml

EXPORT_INCLUDES = $(ML_INCLUDES) $(GALERI_INCLUDES)

AM_CPPFLAGS = $(EXPORT_INCLUDES) $(nox_incdir)

#AM_CPPFLAGS = -I$(srcdir)/Comm \
#	-I$(srcdir)/Operator   \
#	-I$(srcdir)/Smoother   \
#	-I$(srcdir)/Coarsen \
#	-I$(srcdir)/Krylov  \
#	-I$(srcdir)/Main    \
#	-I$(srcdir)/FEGrid  \
#	-I$(srcdir)/Utils   \
#	-I$(srcdir)/MLAPI   \
#	-I$(srcdir)/Include \
#	$(aztecoo_incdir)   \
#	$(ifpack_incdir)    \
#	$(amesos_incdir)    \
#	$(teuchos_incdir)   \
#	$(anasazi_incdir)   \
#	$(triutils_incdir)  \
#	$(epetra_incdir)    \
#        $(epetraext_incdir) \
#        $(nox_incdir)



# ------------------------------------------------------------------------
# ml library specifications
# ------------------------------------------------------------------------
ML_LIB = libml.a

ML_H = \
	$(CORE_H) \
	$(UTIL_H) \
	$(ENUM_H) \
	$(UNIPROC_H) \
	$(SERIAL_H) \
	$(USERBASE_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H) 
#	$(NEWP_SWAHILI_H)

libml_a_SOURCES = \
	$(CORE) \
	$(UTIL) \
	$(UNIPROC) \
	$(SERIAL) \
	$(DIST) \
	$(COMP) \
	$(MPI) 
#	$(NEWP_SWAHILI)

EXTRA_libml_a_SOURCES = \
	$(EXTRA_MPI) \
	$(EXTRA_MPI_H) 
#	$(EXTRA_NEWP_SWAHILI) \
#	$(EXTRA_NEWP_SWAHILI_H)

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libml_a_AR = $(ALTERNATE_AR)
else

libml_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

if ML_USING_TEUCHOS

libml_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

else

libml_a_LIBADD = $(XTRALDADD)

endif

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(ML_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(ML_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: ml_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = ml_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 


# ------------------------------------------------------------------------
# Library Target for Circular Dependencies
# ------------------------------------------------------------------------
if ML_USING_EPETRA
ML_USING_EPETRA_LIBS = cd $(top_builddir)/../epetra ; $(MAKE) all-libs 
else
ML_USING_EPETRA_LIBS = 
endif

if ML_USING_AMESOS
ML_USING_AMESOS_LIBS = cd $(top_builddir)/../amesos ; $(MAKE) all-libs 
else
ML_USING_AMESOS_LIBS = 
endif

if ML_USING_IFPACK
ML_USING_IFPACK_LIBS = cd $(top_builddir)/../ifpack ; $(MAKE) all-libs 
else
ML_USING_IFPACK_LIBS = 
endif

if ML_USING_TEUCHOS
ML_USING_TEUCHOS_LIBS = cd $(top_builddir)/../teuchos ; $(MAKE) all-libs 
else
ML_USING_TEUCHOS_LIBS = 
endif

#if ML_USING_ANASAZI
#ML_USING_ANASAZI_LIBS = cd $(top_builddir)/../anasazi ; $(MAKE) all-libs 
#else
#ML_USING_ANASAZI_LIBS = 
#endif

if ML_USING_AZTECOO
ML_USING_AZTECOO_LIBS = cd $(top_builddir)/../aztecoo ; $(MAKE) all-libs 
else
ML_USING_AZTECOO_LIBS = 
endif

if ML_USING_EPETRAEXT
ML_USING_EPETRAEXT_LIBS = cd $(top_builddir)/../epetraext ; $(MAKE) all-libs 
else
ML_USING_EPETRAEXT_LIBS = 
endif

all-libs:
	$(MAKE) libml.a
	$(ML_USING_EPETRA_LIBS)
	$(ML_USING_AMESOS_LIBS)
	$(ML_USING_IFPACK_LIBS)
	$(ML_USING_TEUCHOS_LIBS)
	$(ML_USING_ANASAZI_LIBS)
	$(ML_USING_AZTECOO_LIBS)
	$(ML_USING_EPETRAEXT_LIBS)
