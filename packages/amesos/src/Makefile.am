# @HEADER
# ************************************************************************
# 
#                 Amesos: Direct Sparse Solver Package
#                 Copyright (2004) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package amesos subdirectory src built successfully."
	@echo ""

# -----------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside Amesos. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

if HAVE_AMESOS_LAPACK
LAPACK_H = $(srcdir)/Amesos_Lapack.h
LAPACK = $(srcdir)/Amesos_Lapack.cpp
else
LAPACK_H =
LAPACK =
endif

if HAVE_AMESOS_SUPERLU
SUPERLU_H = \
	$(srcdir)/Amesos_Superlu.h
SUPERLU = \
	$(srcdir)/Amesos_Superlu.cpp
else
SUPERLU_H =
SUPERLU =
endif

if HAVE_AMESOS_SUPERLUDIST
SUPERLUDIST_H = \
	$(srcdir)/Amesos_Superludist.h
SUPERLUDIST = \
	$(srcdir)/Amesos_Superludist.cpp
else
SUPERLUDIST_H =
SUPERLUDIST =
endif

SUPPORT_H = 
SUPPORT = $(srcdir)/Amesos_Status.cpp  $(srcdir)/Amesos_Control.cpp $(srcdir)/Amesos_Support.cpp


if HAVE_AMESOS_DSCPACK
DSCPACK_H = \
	$(srcdir)/Amesos_Dscpack.h
DSCPACK = \
	$(srcdir)/Amesos_Dscpack.cpp
else
DSCPACK_H =
DSCPACK =
endif

if HAVE_AMESOS_UMFPACK
UMFPACK_H = \
	$(srcdir)/Amesos_Umfpack.h
UMFPACK = \
	$(srcdir)/Amesos_Umfpack.cpp
else
UMFPACK_H =
UMFPACK =
endif

if HAVE_AMESOS_MC64
MC64_H = \
	$(srcdir)/Amesos_MC64.h
MC64 = \
	$(srcdir)/Amesos_MC64.cpp
else
MC64_H =
MC64 =
endif

#if HAVE_AMESOS_SCALAPACK
#SCALAPACK_H = \
#	$(srcdir)/Amesos_Scalapack.h
#SCALAPACK = \
#	$(srcdir)/Amesos_Scalapack.cpp
#else
#SCALAPACK_H =
#SCALAPACK =
#endif

#if HAVE_AMESOS_MERIKOS
#MERIKOS_H = \
#	$(srcdir)/Amesos_Merikos.h
#MERIKOS = \
#	$(srcdir)/Amesos_Merikos.cpp
#else
#MERIKOS_H =
#MERIKOS =
#endif

#-------------------------------------------------------------------------------
#  Here is where we will add Tim Davis' files 
#
if HAVE_AMESOS_KLU
#
#  These files do not need to be copied to the include directory on install,
#  hence they are commented out from here.  
#
KLU_H = \
	$(srcdir)/SuiteSparse/KLU/Include/amesos_klu_decl.h \
	$(srcdir)/SuiteSparse/KLU/Include/amesos_klu_internal.h \
	$(srcdir)/SuiteSparse/KLU/Include/amesos_klu_version.h \
	$(srcdir)/SuiteSparse/BTF/Include/amesos_btf_decl.h \
	$(srcdir)/SuiteSparse/BTF/Include/amesos_btf_internal.h \
	$(srcdir)/Amesos_Klu.h

KLU = \
	$(srcdir)/Amesos_Klu.cpp \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_analyze.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_analyze_given.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_impl.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_defaults.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_diagnostics.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_dump.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_extract.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_factor.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_free_numeric.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_free_symbolic.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_kernel.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_memory.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_refactor.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_scale.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_solve.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_sort.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_tsolve.c \
	$(srcdir)/SuiteSparse/BTF/Source/amesos_btf_order.c \
	$(srcdir)/SuiteSparse/BTF/Source/amesos_btf_maxtrans.c \
	$(srcdir)/SuiteSparse/BTF/Source/amesos_btf_strongcomp.c


else
KLU_H =
KLU =
endif
#-------------------------------------------------------------------------------

if HAVE_AMESOS_MUMPS
MUMPS_H = \
	$(srcdir)/Amesos_Mumps.h
MUMPS = \
	$(srcdir)/Amesos_Mumps.cpp
else
MUMPS_H =
MUMPS =
endif

if HAVE_AMESOS_SCALAPACK
SCALAPACK_H = \
	$(srcdir)/Amesos_Scalapack.h \
	$(srcdir)/Amesos_SCALAPACK_wrappers.h
SCALAPACK = \
	$(srcdir)/Amesos_Scalapack.cpp
else
SCALAPACK_H =
SCALAPACK =
endif

if HAVE_AMESOS_PARDISO
PARDISO_H = \
	$(srcdir)/Amesos_Pardiso.h
PARDISO = \
	$(srcdir)/Amesos_Pardiso.cpp
else
PARDISO_H =
PARDISO =
endif

if HAVE_AMESOS_TAUCS
TAUCS_H = \
	$(srcdir)/Amesos_Taucs.h
TAUCS = \
	$(srcdir)/Amesos_Taucs.cpp
else
TAUCS_H =
TAUCS =
endif


PARAKLETE_KLU_COMMON = \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_1.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_2.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_aat.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_control.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_defaults.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_dump.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_global.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_info.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_order.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_postorder.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_post_tree.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_preprocess.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_valid.c \
	$(srcdir)/SuiteSparse/COLAMD/Source/amesos_colamd.c \
	$(srcdir)/SuiteSparse/COLAMD/Source/amesos_colamd_global.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_common.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_complex.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_error.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_memory.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_sparse.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_transpose.c \
	$(srcdir)/SuiteSparse/CAMD/Source/amesos_camd_2.c \
	$(srcdir)/SuiteSparse/CAMD/Source/amesos_camd_postorder.c

PARAKLETE_KLU_COMMON_H = \
	$(srcdir)/SuiteSparse/AMD/Include/amesos_amd.h \
	$(srcdir)/SuiteSparse/AMD/Include/amesos_amd_internal.h \
	$(srcdir)/SuiteSparse/COLAMD/Include/amesos_colamd.h \
	$(srcdir)/SuiteSparse/CHOLMOD/Include/amesos_cholmod_check.h \
	$(srcdir)/SuiteSparse/CHOLMOD/Include/amesos_cholmod_core.h \
	$(srcdir)/SuiteSparse/CHOLMOD/Include/amesos_cholmod.h \
	$(srcdir)/SuiteSparse/CHOLMOD/Include/amesos_cholmod_internal.h \
	$(srcdir)/SuiteSparse/CHOLMOD/Include/amesos_cholmod_io64.h \
	$(srcdir)/SuiteSparse/CHOLMOD/Include/amesos_cholmod_partition.h \
	$(srcdir)/SuiteSparse/UFconfig/amesos_UFconfig.h

if HAVE_AMESOS_PARAKLETE
PARAKLETE_KLU = ${PARAKLETE_KLU_COMMON}
PARAKLETE_KLU_H = ${PARAKLETE_KLU_COMMON_H}
else
if HAVE_AMESOS_KLU
PARAKLETE_KLU = ${PARAKLETE_KLU_COMMON}
PARAKLETE_KLU_H = ${PARAKLETE_KLU_COMMON_H}
else
PARAKLETE_KLU = 
endif 
endif

if HAVE_AMESOS_PARAKLETE
PARAKLETE_H = \
	$(srcdir)/Amesos_Paraklete.h \
	$(srcdir)/SuiteSparse/PARAKLETE/Include/amesos_paraklete_decl.h \
	$(srcdir)/SuiteSparse/CHOLMOD/Include/amesos_cholmod_complexity.h \
	$(srcdir)/SuiteSparse/CCOLAMD/Include/amesos_ccolamd.h \
	$(srcdir)/SuiteSparse/CHOLMOD/Include/amesos_cholmod_blas.h \
	$(srcdir)/SuiteSparse/CHOLMOD/Include/amesos_cholmod_cholesky.h \
	$(srcdir)/SuiteSparse/CHOLMOD/Include/amesos_cholmod_config.h \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_t_cholmod_dense.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Include/amesos_cholmod_template.h \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_t_cholmod_change_factor.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_t_cholmod_transpose.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_t_cholmod_triplet.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Cholesky/amesos_t_cholmod_solve.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Cholesky/amesos_t_cholmod_lsolve.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Cholesky/amesos_t_cholmod_ltsolve.c \
	$(srcdir)/SuiteSparse/CAMD/Include/amesos_camd.h \
	$(srcdir)/SuiteSparse/CAMD/Include/amesos_camd_internal.h 

#  Note:  PARAKLETE requires KLU (for colamd) 

PARAKLETE = \
	$(srcdir)/SuiteSparse/PARAKLETE/Source/amesos_paraklete_analyze.c \
	$(srcdir)/SuiteSparse/PARAKLETE/Source/amesos_paraklete_btf.c \
	$(srcdir)/SuiteSparse/PARAKLETE/Source/amesos_paraklete_error.c \
	$(srcdir)/SuiteSparse/PARAKLETE/Source/amesos_paraklete_factorize.c \
	$(srcdir)/SuiteSparse/PARAKLETE/Source/amesos_paraklete_finish.c \
	$(srcdir)/SuiteSparse/PARAKLETE/Source/amesos_paraklete_kernel.c \
	$(srcdir)/SuiteSparse/PARAKLETE/Source/amesos_paraklete_lsolve_node.c \
	$(srcdir)/SuiteSparse/PARAKLETE/Source/amesos_paraklete_node.c \
	$(srcdir)/SuiteSparse/PARAKLETE/Source/amesos_paraklete_reanalyze.c \
	$(srcdir)/SuiteSparse/PARAKLETE/Source/amesos_paraklete_solve.c \
	$(srcdir)/SuiteSparse/PARAKLETE/Source/amesos_paraklete_start.c \
	$(srcdir)/SuiteSparse/PARAKLETE/Source/amesos_paraklete_usolve_node.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_common.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_dense.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_factor.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_change_factor.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_memory.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_complex.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_transpose.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_sparse.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_band.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_copy.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_triplet.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_error.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_aat.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_add.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Check/amesos_cholmod_check.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_l_postorder.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_l_amd.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_l_colamd.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_l_etree.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_l_rowcolcounts.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_l_analyze.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Partition/amesos_cholmod_l_ccolamd.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Partition/amesos_cholmod_l_camd.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Partition/amesos_cholmod_l_csymamd.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Partition/amesos_cholmod_l_metis.c \
	$(srcdir)/SuiteSparse/CHOLMOD/Partition/amesos_cholmod_l_nesdis.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_l1.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_l2.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_l_aat.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_l_order.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_l_postorder.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_l_post_tree.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_l_preprocess.c \
	$(srcdir)/SuiteSparse/AMD/Source/amesos_amd_l_valid.c \
	$(srcdir)/SuiteSparse/CAMD/Source/amesos_camd_global.c \
	$(srcdir)/SuiteSparse/CAMD/Source/amesos_camd_l2.c \
	$(srcdir)/SuiteSparse/CAMD/Source/amesos_camd_l_postorder.c \
	$(srcdir)/SuiteSparse/COLAMD/Source/amesos_colamd_l.c \
	$(srcdir)/SuiteSparse/CCOLAMD/Source/amesos_ccolamd_l.c \
	$(srcdir)/SuiteSparse/CCOLAMD/Source/amesos_ccolamd_global.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_l.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_l_analyze.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_l_analyze_given.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_l_defaults.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_l_free_numeric.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_l_free_symbolic.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_l_factor.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_l_kernel.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_l_memory.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_l_scale.c \
	$(srcdir)/SuiteSparse/KLU/Source/amesos_klu_l_solve.c \
	$(srcdir)/SuiteSparse/BTF/Source/amesos_btf_l_order.c \
	$(srcdir)/SuiteSparse/BTF/Source/amesos_btf_l_maxtrans.c \
	$(srcdir)/SuiteSparse/BTF/Source/amesos_btf_l_strongcomp.c \
	$(srcdir)/Amesos_Paraklete.cpp
else
PARAKLETE_H =
PARAKLETE =
endif

#
#  This file does not need to be copied to the include directory
#	$(srcdir)/CrsMatrixTranspose.h
COMP_H = \
	$(srcdir)/Amesos_BaseSolver.h \
	$(srcdir)/Amesos_ComponentBaseSolver.h \
	$(srcdir)/Amesos_ConfigDefs.h \
	$(srcdir)/CrsMatrixTranspose.h  \
        $(srcdir)/Amesos_NoCopiable.h \
        $(srcdir)/Amesos_Utils.h \
        $(srcdir)/Amesos_Time.h \
        $(srcdir)/Amesos_Control.h \
        $(srcdir)/Amesos_Support.h \
        $(srcdir)/Amesos_Status.h \
	$(srcdir)/Amesos_Scaling.h \
	$(srcdir)/Amesos_Reordering.h

COMP = \
	$(srcdir)/CrsMatrixTranspose.cpp

FACTORY_H = \
	$(srcdir)/Amesos.h
FACTORY = \
	$(srcdir)/Amesos.cpp

#
#  I turned off the addition of external/slu to the include path 
#  so that I could get SLU dist to work.  Ken Feb 11, 2003
#
#Specify locations of include files within Trilinos, but outside Amesos
# AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src \
#		-I$(top_srcdir)/../teuchos/src -I$(top_builddir)/../teuchos/src \
#	       -I$(top_srcdir)/../../external/slu -I$(top_srcdir)/src/ksparse

# if BUILD_SLUSERIAL
# SLUSERIAL_CPPFLAGS = -I$(top_srcdir)/../../external/slu
# else
# SLUSERIAL_CPPFLAGS = 
# endif

include $(top_builddir)/Makefile.export.amesos

if HAVE_AMESOS_EPETRAEXT
   EPETRA_EXT_INCLUDE =	-I$(top_srcdir)/../epetraext/src/transform -I$(top_builddir)/../epetraext/src/transform -I$(top_srcdir)/../epetraext/src -I$(top_builddir)/../epetraext/src 
else
   EPETRA_EXT_INCLUDE = 
endif 

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(AMESOS_INCLUDES))
else
EXPORT_INCLUDES = $(AMESOS_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES) $(EPETRA_EXT_INCLUDE) $(SLUSERIAL_CPPFLAGS) \
	      -I$(top_srcdir)/src/SuiteSparse/AMD/Include \
	      -I$(top_srcdir)/src/SuiteSparse/PARAKLETE/Include \
	      -I$(top_srcdir)/src/SuiteSparse/CHOLMOD/Include \
	      -I$(top_srcdir)/src/SuiteSparse/CCOLAMD/Include \
	      -I$(top_srcdir)/src/SuiteSparse/COLAMD/Include \
	      -I$(top_srcdir)/src/SuiteSparse/BTF/Include \
	      -I$(top_srcdir)/src/SuiteSparse/CAMD/Include \
	      -I$(top_srcdir)/src/SuiteSparse/UFconfig \
	      -I$(top_srcdir)/src/SuiteSparse/KLU/Include 

# ------------------------------------------------------------------------
# Amesos library specifications
# ------------------------------------------------------------------------

AMESOS_LIB = libamesos.a


AMESOS_H = \
	$(SUPPORT_H) \
	$(DSCPACK_H) \
	$(UMFPACK_H) \
	$(MUMPS_H) \
	$(PARDISO_H) \
	$(TAUCS_H) \
	$(PARAKLETE_KLU_H) \
	$(KLU_H) \
	$(PARAKLETE_H) \
	$(MC64_H) \
	$(SCALAPACK_H) \
	$(FACTORY_H) \
	$(SUPERLU_H) \
	$(SUPERLUDIST_H) \
	$(LAPACK_H) \
	$(COMP_H) \
	$(srcdir)/Amesos_BaseSolver.h \
	$(srcdir)/Amesos_ComponentBaseSolver.h \
	$(srcdir)/Amesos_TestRowMatrix.h


libamesos_a_SOURCES = \
	$(SUPPORT) \
	$(DSCPACK) \
	$(UMFPACK) \
	$(MUMPS) \
	$(PARDISO) \
	$(TAUCS) \
	$(PARAKLETE_KLU) \
	$(KLU) \
	$(PARAKLETE) \
	$(MC64) \
	$(SCALAPACK) \
	$(FACTORY) \
	$(SUPERLU) \
	$(SUPERLUDIST) \
	$(LAPACK) \
	$(COMP)

EXTRA_libamesos_a_SOURCES = \
	$(EXTRA_MPI) \
	$(EXTRA_MPI_H) \
	$(EXTRA_MPISMP) \
	$(EXTRA_MPISMP_H)

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libamesos_a_AR = $(ALTERNATE_AR)
else

libamesos_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libamesos_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(AMESOS_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(AMESOS_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Amesos_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(includedir)
nodist_execinclude_HEADERS = Amesos_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 
