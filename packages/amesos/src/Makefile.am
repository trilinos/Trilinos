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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
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
	$(srcdir)/UFsparse/KLU/Include/klu.h \
	$(srcdir)/UFsparse/KLU/Include/klu_internal.h \
	$(srcdir)/UFsparse/KLU/Include/klu_version.h \
	$(srcdir)/UFsparse/BTF/Include/btf.h \
	$(srcdir)/UFsparse/BTF/Include/btf_internal.h \
	$(srcdir)/Amesos_Klu.h

KLU = \
	$(srcdir)/Amesos_Klu.cpp \
	$(srcdir)/UFsparse/KLU/Source/klu_analyze.c \
	$(srcdir)/UFsparse/KLU/Source/klu_analyze_given.c \
	$(srcdir)/UFsparse/KLU/Source/klu.c \
	$(srcdir)/UFsparse/KLU/Source/klu_defaults.c \
	$(srcdir)/UFsparse/KLU/Source/klu_diagnostics.c \
	$(srcdir)/UFsparse/KLU/Source/klu_dump.c \
	$(srcdir)/UFsparse/KLU/Source/klu_extract.c \
	$(srcdir)/UFsparse/KLU/Source/klu_factor.c \
	$(srcdir)/UFsparse/KLU/Source/klu_free_numeric.c \
	$(srcdir)/UFsparse/KLU/Source/klu_free_symbolic.c \
	$(srcdir)/UFsparse/KLU/Source/klu_kernel.c \
	$(srcdir)/UFsparse/KLU/Source/klu_memory.c \
	$(srcdir)/UFsparse/KLU/Source/klu_refactor.c \
	$(srcdir)/UFsparse/KLU/Source/klu_scale.c \
	$(srcdir)/UFsparse/KLU/Source/klu_solve.c \
	$(srcdir)/UFsparse/KLU/Source/klu_sort.c \
	$(srcdir)/UFsparse/KLU/Source/klu_tsolve.c \
	$(srcdir)/UFsparse/BTF/Source/btf_order.c \
	$(srcdir)/UFsparse/BTF/Source/maxtrans.c \
	$(srcdir)/UFsparse/BTF/Source/strongcomp.c


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

if HAVE_AMESOS_PASTIX
PASTIX_H = \
	$(srcdir)/Amesos_Pastix.h
PASTIX = \
	$(srcdir)/Amesos_Pastix.cpp
else
PASTIX_H =
PASTIX =
endif


PARAKLETE_KLU_COMMON = \
	$(srcdir)/UFsparse/AMD/Source/amd_1.c \
	$(srcdir)/UFsparse/AMD/Source/amd_2.c \
	$(srcdir)/UFsparse/AMD/Source/amd_aat.c \
	$(srcdir)/UFsparse/AMD/Source/amd_control.c \
	$(srcdir)/UFsparse/AMD/Source/amd_defaults.c \
	$(srcdir)/UFsparse/AMD/Source/amd_dump.c \
	$(srcdir)/UFsparse/AMD/Source/amd_global.c \
	$(srcdir)/UFsparse/AMD/Source/amd_info.c \
	$(srcdir)/UFsparse/AMD/Source/amd_order.c \
	$(srcdir)/UFsparse/AMD/Source/amd_postorder.c \
	$(srcdir)/UFsparse/AMD/Source/amd_post_tree.c \
	$(srcdir)/UFsparse/AMD/Source/amd_preprocess.c \
	$(srcdir)/UFsparse/AMD/Source/amd_valid.c \
	$(srcdir)/UFsparse/COLAMD/colamd.c \
	$(srcdir)/UFsparse/COLAMD/colamd_global.c 

PARAKLETE_KLU_COMMON_H = \
	$(srcdir)/UFsparse/AMD/Include/amd.h \
	$(srcdir)/UFsparse/AMD/Include/amd_internal.h \
	$(srcdir)/UFsparse/COLAMD/colamd.h \
	$(srcdir)/UFsparse/CHOLMOD/Include/cholmod_check.h \
	$(srcdir)/UFsparse/CHOLMOD/Include/cholmod_core.h \
	$(srcdir)/UFsparse/CHOLMOD/Include/cholmod.h \
	$(srcdir)/UFsparse/CHOLMOD/Include/cholmod_internal.h \
	$(srcdir)/UFsparse/CHOLMOD/Include/cholmod_partition.h \
	$(srcdir)/UFsparse/UFconfig/UFconfig.h

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
	$(srcdir)/UFsparse/PARAKLETE/paraklete.h \
	$(srcdir)/UFsparse/CCOLAMD/ccolamd.h 

#  Note:  PARAKLETE rewuires KLU (for colamd) 

BOGUSMEX= \
	$(srcdir)/UFsparse/CCOLAMD/ccolamdmex.c \
	$(srcdir)/UFsparse/CCOLAMD/csymamdmex.c

PARAKLETE = \
	$(srcdir)/UFsparse/PARAKLETE/memory.c \
	$(srcdir)/UFsparse/PARAKLETE/paraklete_analyze.c \
	$(srcdir)/UFsparse/PARAKLETE/paraklete_factorize.c \
	$(srcdir)/UFsparse/PARAKLETE/paraklete_kernel.c \
	$(srcdir)/UFsparse/PARAKLETE/paraklete_solve.c \
	$(srcdir)/UFsparse/PARAKLETE/paraklete_node.c \
	$(srcdir)/UFsparse/PARAKLETE/paraklete_lsolve_node.c \
	$(srcdir)/UFsparse/PARAKLETE/paraklete_usolve_node.c \
	$(srcdir)/UFsparse/CHOLMOD/Core/cholmod_common.c \
	$(srcdir)/UFsparse/CHOLMOD/Core/cholmod_dense.c \
	$(srcdir)/UFsparse/CHOLMOD/Core/cholmod_factor.c \
	$(srcdir)/UFsparse/CHOLMOD/Core/cholmod_change_factor.c \
	$(srcdir)/UFsparse/CHOLMOD/Core/cholmod_memory.c \
	$(srcdir)/UFsparse/CHOLMOD/Core/cholmod_complex.c \
	$(srcdir)/UFsparse/CHOLMOD/Core/cholmod_transpose.c \
	$(srcdir)/UFsparse/CHOLMOD/Core/cholmod_sparse.c \
	$(srcdir)/UFsparse/CHOLMOD/Core/cholmod_band.c \
	$(srcdir)/UFsparse/CHOLMOD/Core/cholmod_copy.c \
	$(srcdir)/UFsparse/CHOLMOD/Core/cholmod_triplet.c \
	$(srcdir)/UFsparse/CHOLMOD/Core/cholmod_error.c \
	$(srcdir)/UFsparse/CHOLMOD/Core/cholmod_aat.c \
	$(srcdir)/UFsparse/CHOLMOD/Core/cholmod_add.c \
	$(srcdir)/UFsparse/CHOLMOD/Check/cholmod_check.c \
	$(srcdir)/UFsparse/CHOLMOD/Cholesky/cholmod_postorder.c \
	$(srcdir)/UFsparse/CHOLMOD/Cholesky/cholmod_amd.c \
	$(srcdir)/UFsparse/CHOLMOD/Cholesky/cholmod_colamd.c \
	$(srcdir)/UFsparse/CHOLMOD/Cholesky/cholmod_etree.c \
	$(srcdir)/UFsparse/CHOLMOD/Cholesky/cholmod_rowcolcounts.c \
	$(srcdir)/UFsparse/CHOLMOD/Cholesky/cholmod_solve.c \
	$(srcdir)/UFsparse/CHOLMOD/Cholesky/cholmod_resymbol.c \
	$(srcdir)/UFsparse/CHOLMOD/Cholesky/cholmod_analyze.c \
	$(srcdir)/UFsparse/CHOLMOD/Partition/cholmod_ccolamd.c \
	$(srcdir)/UFsparse/CHOLMOD/Partition/cholmod_camd.c \
	$(srcdir)/UFsparse/CHOLMOD/Partition/cholmod_csymamd.c \
	$(srcdir)/UFsparse/CHOLMOD/Partition/cholmod_metis.c \
	$(srcdir)/UFsparse/CHOLMOD/Partition/cholmod_nesdis.c \
	$(srcdir)/UFsparse/CAMD/Source/camd_1.c \
	$(srcdir)/UFsparse/CAMD/Source/camd_2.c \
	$(srcdir)/UFsparse/CAMD/Source/camd_aat.c \
	$(srcdir)/UFsparse/CAMD/Source/camd_control.c \
	$(srcdir)/UFsparse/CAMD/Source/camd_defaults.c \
	$(srcdir)/UFsparse/CAMD/Source/camd_dump.c \
	$(srcdir)/UFsparse/CAMD/Source/camd_info.c \
	$(srcdir)/UFsparse/CAMD/Source/camd_order.c \
	$(srcdir)/UFsparse/CAMD/Source/camd_postorder.c \
	$(srcdir)/UFsparse/CAMD/Source/camd_preprocess.c \
	$(srcdir)/UFsparse/CAMD/Source/camd_valid.c \
	$(srcdir)/UFsparse/CAMD/Source/camd_global.c \
	$(srcdir)/UFsparse/CCOLAMD/ccolamd.c \
	$(srcdir)/UFsparse/CCOLAMD/ccolamd_global.c \
	$(srcdir)/Amesos_Paraklete.cpp
else
PARAKLETE_H =
PARAKLETE =
endif




if HAVE_AMESOS_KUNDERT
KUNDERT_H = \
       $(srcdir)/Epetra_CrsKundertSparse.h
KUNDERT = \
       $(srcdir)/Epetra_CrsKundertSparse.cpp \
       $(srcdir)/ksparse/alloc.c \
       $(srcdir)/ksparse/smplinkrows.c \
       $(srcdir)/ksparse/spbuild.c \
       $(srcdir)/ksparse/spoutput.c \
       $(srcdir)/ksparse/spsolve.c \
       $(srcdir)/ksparse/spalloc.c \
       $(srcdir)/ksparse/spfactor.c \
       $(srcdir)/ksparse/spsmp.c \
       $(srcdir)/ksparse/spmatrix.h \
       $(srcdir)/ksparse/sperror.h \
       $(srcdir)/ksparse/spconfig.h \
       $(srcdir)/ksparse/iferrmsg.h \
       $(srcdir)/ksparse/strext.h \
       $(srcdir)/ksparse/spice.h \
       $(srcdir)/ksparse/util.h \
       $(srcdir)/ksparse/misc.h \
       $(srcdir)/ksparse/suffix.h \
       $(srcdir)/ksparse/port.h \
       $(srcdir)/ksparse/hw.h \
       $(srcdir)/ksparse/config.h \
       $(srcdir)/ksparse/capabil.h \
       $(srcdir)/ksparse/os_linux.h \
       $(srcdir)/ksparse/os_unix.h \
       $(srcdir)/ksparse/hw_ieee.h \
       $(srcdir)/ksparse/smpdefs.h \
       $(srcdir)/ksparse/spdefs.h \
       $(srcdir)/ksparse/sppars.h \
       $(srcdir)/ksparse/spmat.h \
       $(srcdir)/ksparse/smpdefs.h \
       $(srcdir)/ksparse/complex.h \
       $(srcdir)/ksparse/sputils.c
else
KUNDERT_H = 
KUNDERT =
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

if HAVE_AMESOS_KUNDERT
KUNDERT_CPPFLAGS = -I$(top_srcdir)/src/ksparse
else
KUNDERT_CPPFLAGS = 
endif

#AM_CPPFLAGS = -I$(top_srcdir)/../epetraext/src/transform \
#	-I$(top_srcdir)/../epetraext/src \
#	-I$(top_srcdir)/../epetra/src \
#	-I$(top_srcdir)/src/UFsparse/AMD/Include \
#	-I$(top_srcdir)/src/UFsparse/AMD/Source \
#	-I$(top_srcdir)/src/colamd \
#	-I$(top_srcdir)/src/btf \
#	-I$(top_srcdir)/src/klu_btf \
#	-I$(top_srcdir)/src/klu_kernel \
#	-I$(top_srcdir)/../teuchos/src  -I$(top_builddir)/../teuchos/src \
#	-I$(top_builddir)/../epetraext/src \
#	-I$(top_builddir)/../epetra/src \
#	      $(SLUSERIAL_CPPFLAGS) $(KUNDERT_CPPFLAGS)

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

UFSPARSE_CPPFLAGS = -DNSUPERNODAL

AM_CPPFLAGS = $(EXPORT_INCLUDES) $(EPETRA_EXT_INCLUDE) $(SLUSERIAL_CPPFLAGS) $(KUNDERT_CPPFLAGS) $(UFSPARSE_CPPFLAGS) \
	      -I$(top_srcdir)/src/UFsparse/AMD/Include \
	      -I$(top_srcdir)/src/UFsparse/PARAKLETE \
	      -I$(top_srcdir)/src/UFsparse/CHOLMOD/Include \
	      -I$(top_srcdir)/src/UFsparse/CCOLAMD \
	      -I$(top_srcdir)/src/UFsparse/COLAMD \
	      -I$(top_srcdir)/src/UFsparse/BTF/Include \
	      -I$(top_srcdir)/src/UFsparse/CAMD/Include \
	      -I$(top_srcdir)/src/UFsparse/UFconfig \
	      -I$(top_srcdir)/src/UFsparse/KLU/Include 

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
	$(PASTIX_H) \
	$(PARAKLETE_KLU_H) \
	$(KLU_H) \
	$(PARAKLETE_H) \
	$(MC64_H) \
	$(SCALAPACK_H) \
	$(FACTORY_H) \
	$(SUPERLU_H) \
	$(SUPERLUDIST_H) \
	$(KUNDERT_H) \
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
	$(PASTIX) \
	$(PARAKLETE_KLU) \
	$(KLU) \
	$(PARAKLETE) \
	$(MC64) \
	$(SCALAPACK) \
	$(FACTORY) \
	$(SUPERLU) \
	$(SUPERLUDIST) \
	$(LAPACK) \
	$(KUNDERT) \
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

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = Amesos_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Library Target for Circular Dependencies
# ------------------------------------------------------------------------
all-libs:
	$(MAKE) libamesos.a
	cd $(top_builddir)/../teuchos ; $(MAKE) all-libs 
	cd $(top_builddir)/../epetra ; $(MAKE) all-libs 
if HAVE_AMESOS_TRIUTILS
	cd $(top_builddir)/../triutils ; $(MAKE) all-libs 
endif
if HAVE_AMESOS_EPETRAEXT
	cd $(top_builddir)/../epetraext ; $(MAKE) all-libs 
endif





