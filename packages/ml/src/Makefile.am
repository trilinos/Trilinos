# @HEADER
# @HEADER

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
$(srcdir)/Coarsen/ml_check.h \
$(srcdir)/Coarsen/ml_aggregate.h \
$(srcdir)/Coarsen/ml_ggraph.h \
$(srcdir)/FEGrid/ml_mapper.h \
$(srcdir)/FEGrid/ml_gridagx.h \
$(srcdir)/FEGrid/ml_pde.h       $(srcdir)/FEGrid/ml_bdrypts.h \
$(srcdir)/Coarsen/ml_agg_reitzinger.h \
$(srcdir)/Coarsen/ml_agg_METIS.h \
$(srcdir)/Coarsen/ml_agg_ParMETIS.h \
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
$(srcdir)/Coarsen/ml_amg_MIS.c \
$(srcdir)/Coarsen/ml_agg_MIS.c \
$(srcdir)/Coarsen/ml_agg_METIS.c \
$(srcdir)/Coarsen/ml_agg_ParMETIS.c \
$(srcdir)/Coarsen/ml_aggregate.c \
$(srcdir)/Coarsen/ml_ggraph.c \
$(srcdir)/Coarsen/ml_agg_reitzinger.c  \
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
$(srcdir)/Utils/ml_aztec_utils.h $(srcdir)/Utils/ml_epetra_utils.h \
$(srcdir)/Utils/ml_utils.h \
$(srcdir)/Utils/ml_memory.h \
$(srcdir)/Utils/ml_agg_info.h \
$(srcdir)/Utils/ml_aztec_lapack.h      $(srcdir)/Utils/ml_rbm.h \
$(srcdir)/Utils/ml_read_utils.h \
$(srcdir)/Utils/ml_superlu_lapack.h \
$(srcdir)/Utils/ml_epetra_operator.h \
$(srcdir)/Utils/ml_twogrid_analysis.h \
$(srcdir)/Utils/ml_intlist.h \
$(srcdir)/Utils/ml_vec.h \
$(srcdir)/Utils/ml_lapack.h            $(srcdir)/Utils/ml_vendor_lapack.h \
$(srcdir)/Utils/ml_utils.h \
$(srcdir)/Utils/ml_eigf2c.h \
$(srcdir)/Utils/ml_vampir.h \
$(srcdir)/Utils/ml_viz_opendx.h \
$(srcdir)/Utils/ml_amesos.h \
$(srcdir)/Utils/ml_amesos_wrap.h \
$(srcdir)/Utils/ml_ifpack.h \
$(srcdir)/Utils/ml_ifpack_wrap.h \
$(srcdir)/Utils/ml_epetra_preconditioner.h \
$(srcdir)/Utils/ml_anasazi.h


UTIL = \
$(srcdir)/Utils/ml_rbm.c \
$(srcdir)/Utils/ml_aztec_utils.c       $(srcdir)/Utils/ml_read_utils.c \
$(srcdir)/Utils/ml_epetra_operator.cpp \
$(srcdir)/Utils/ml_twogrid_analysis.c \
$(srcdir)/Utils/ml_epetra_utils.cpp \
$(srcdir)/Utils/ml_utils.c             $(srcdir)/Utils/ml_intlist.c \
$(srcdir)/Utils/ml_vec.c \
$(srcdir)/Utils/ml_lapack.c \
$(srcdir)/Utils/ml_memory.c            $(srcdir)/Utils/ml_eigf2c.c \
$(srcdir)/Utils/ml_agg_info.c \
$(srcdir)/Utils/ml_viz_opendx.c \
$(srcdir)/Utils/ml_amesos.c \
$(srcdir)/Utils/ml_amesos_wrap.cpp \
$(srcdir)/Utils/ml_ifpack.c \
$(srcdir)/Utils/ml_ifpack_wrap.cpp \
$(srcdir)/Utils/ml_MultiLevelPreconditioner.cpp \
$(srcdir)/Utils/ml_MultiLevelPreconditioner_Filtering.cpp \
$(srcdir)/Utils/ml_MultiLevelPreconditioner_NullSpace.cpp \
$(srcdir)/Utils/ml_MultiLevelPreconditioner_Defaults.cpp \
$(srcdir)/Utils/ml_anasazi.cpp


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
#  
AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src \
              -I$(top_srcdir)/../aztecoo/src -I$(top_builddir)/../aztecoo/src \
              -I$(top_srcdir)/../amesos/src -I$(top_builddir)/../amesos/src \
              -I$(top_srcdir)/../ifpack/src -I$(top_builddir)/../ifpack/src \
              -I$(top_srcdir)/../teuchos/src -I$(top_builddir)/../teuchos/src \
              -I$(top_srcdir)/../anasazi/src -I$(top_builddir)/../anasazi/src \
              -I$(top_srcdir)/../triutils/src -I$(top_builddir)/../triutils/src \
              -I$(srcdir)/Comm \
              -I$(srcdir)/Operator \
              -I$(srcdir)/Smoother \
              -I$(srcdir)/Coarsen \
              -I$(srcdir)/Krylov \
              -I$(srcdir)/Main \
              -I$(srcdir)/FEGrid \
              -I$(srcdir)/Utils \
              -I$(srcdir)/Include



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






