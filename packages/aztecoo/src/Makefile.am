# @HEADER
# ************************************************************************
# 
#        AztecOO: An Object-Oriented Aztec Linear Solver Package 
#                 Copyright (2002) Sandia Corporation
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

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside AztecOO. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#WITH_CPP category is not incorporated into the make structure yet.  Must
#make an option to use it or not with default being yes.  Must also 
#incorporate platform specific files after I figure out exactly how to 
# define a platform for Aztec.  Also change library name if not compiling
# c++ code.

# I don't know what to do with the following files - many are platform
# specific - note also that some of the files not listed in the
# following group may also be specific to the Linux platform
# the group includes the files listed under mpi

#	$(srcdir)/az_comm_ncube.c 
#	$(srcdir)/az_comm_intel.c 

CORE_H = \
	$(srcdir)/az_aztec_defs.h \
	$(srcdir)/az_aztecf.h \
	$(srcdir)/az_aztec.h \
	$(srcdir)/az_f77func.h \
	$(srcdir)/az_y12m_wrappers.h \
	$(srcdir)/az_blas_wrappers.h \
	$(srcdir)/az_lapack_wrappers.h

CORE = \
	$(srcdir)/az_comm.c 
# The above file might be platform specific

### Eventually make the following an option (default is to include)
### Making this optional is one step toward removing the current
### dependency on Epetra when only Aztec is needed
AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src \
	-I$(top_srcdir)/../ifpack/src -I$(top_builddir)../ifpack/src
WITH_CPP_H = \
	$(srcdir)/AztecOO.h \
	$(srcdir)/AztecOO_Operator.h \
	$(srcdir)/Epetra_MsrMatrix.h \
	$(srcdir)/Aztec2Petra.h \
	$(srcdir)/AZOO_iterate.h \
        $(srcdir)/AztecOO_StatusTest.h \
        $(srcdir)/AztecOO_StatusType.h \
        $(srcdir)/AztecOO_StatusTestCombo.h \
        $(srcdir)/AztecOO_StatusTestMaxIters.h \
        $(srcdir)/AztecOO_StatusTestResNorm.h


WITH_CPP = \
	$(srcdir)/AztecOO.cpp \
	$(srcdir)/AztecOO_Operator.cpp \
	$(srcdir)/Epetra_MsrMatrix.cpp \
	$(srcdir)/Aztec2Petra.cpp \
	$(srcdir)/AZOO_iterate.cpp \
        $(srcdir)/AztecOO_StatusTestCombo.cpp \
        $(srcdir)/AztecOO_StatusTestMaxIters.cpp \
        $(srcdir)/AztecOO_StatusTestResNorm.cpp

### end option 
UTIL_H = \
	$(srcdir)/AztecOO_Version.h 

UTIL = \
	$(srcdir)/az_check.c \
	$(srcdir)/az_flop_cnt.c \
	$(srcdir)/az_ilu_util.c \
	$(srcdir)/az_interface.c \
	$(srcdir)/az_matrix_util.c \
	$(srcdir)/az_pad_utils.c \
	$(srcdir)/az_tools.c \
	$(srcdir)/az_f_util.f \
	$(srcdir)/az_util.c

ENUM_H = 

UNIPROC_H = 

UNIPROC = 

EXTRA_SERIAL_H = 

EXTRA_SERIAL = \
$(srcdir)/md_timer_generic.c \
	$(srcdir)/md_wrap_scalar_c.c 


EXTRA_MPI_H = 

EXTRA_MPI = \
	$(srcdir)/md_timer_mpi.c \
	$(srcdir)/md_wrap_mpi_c.c 

EXTRA_DIST = \
	$(EXTRA_MPI_H) \
	$(EXTRA_MPI) \
	$(EXTRA_SERIAL_H) \
	$(EXTRA_SERIAL)

if HAVE_MPI

MPI_H = \
	$(EXTRA_MPI_H)
MPI = \
	$(EXTRA_MPI)

SERIAL_H = 
SERIAL = 

else

MPI_H =
MPI =

SERIAL_H = \
	$(EXTRA_SERIAL_H)

SERIAL = \
	$(EXTRA_SERIAL)

endif

#Aztec does not yet support threads
#EXTRA_MPISMP_H =

#EXTRA_MPISMP =

#if HAVE_THREADS

#MPISMP_H = \
#	$(EXTRA_MPISMP_H)
#MPISMP = \
#	$(EXTRA_MPISMP)

#else

MPISMP_H =
MPISMP =

#endif

USERBASE_H = 


LANG_WRAPS_H = 

LANG_WRAPS = \
	$(srcdir)/az_fortran_wrap.c

DIST_H = 

DIST = \
	$(srcdir)/az_bilu.c \
	$(srcdir)/az_cg.c \
	$(srcdir)/az_cgs.c \
	$(srcdir)/az_cgstab.c \
	$(srcdir)/az_converge.c \
	$(srcdir)/az_dd_overlap.c \
	$(srcdir)/az_dgemv2.c \
	$(srcdir)/az_dgemv3.c \
	$(srcdir)/az_domain_decomp.c \
	$(srcdir)/az_fix_pt.c \
	$(srcdir)/az_gmres.c \
	$(srcdir)/az_gmresr.c \
	$(srcdir)/az_gsumd_puma.c \
	$(srcdir)/az_icc.c \
	$(srcdir)/az_ilut.c \
	$(srcdir)/az_lu_y12.c \
	$(srcdir)/az_poly.c \
	$(srcdir)/az_qmrcgs.c \
	$(srcdir)/az_rilu.c \
	$(srcdir)/az_subdomain_solver.c \
	$(srcdir)/az_cg_condnum.c \
	$(srcdir)/az_gmres_condnum.c

#There is no evidence that this file is used by anyone
#	$(srcdir)/az_xyce_gmres.c

COMP_H = 

COMP = \
	$(srcdir)/az_colperm.c \
	$(srcdir)/az_matvec_mult.c \
	$(srcdir)/az_old_matvec_mult.c \
	$(srcdir)/az_precond.c \
	$(srcdir)/az_reorder.f \
	$(srcdir)/az_scaling.c \
	$(srcdir)/az_solve.c \
	$(srcdir)/az_sort.c 


# ------------------------------------------------------------------------
# Aztec library specifications
# ------------------------------------------------------------------------

AZTECOO_LIB = libaztecoo.a


AZTECOO_H = \
	$(CORE_H) \
	$(WITH_CPP_H) \
	$(UTIL_H) \
	$(ENUM_H) \
	$(UNIPROC_H) \
	$(SERIAL_H) \
	$(USERBASE_H) \
	$(LANG_WRAPS_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H) \
	$(MPISMP_H)

libaztecoo_a_SOURCES = \
	$(CORE) \
	$(WITH_CPP) \
	$(UTIL) \
	$(UNIPROC) \
	$(SERIAL) \
	$(LANG_WRAPS) \
	$(DIST) \
	$(COMP) \
	$(MPI) \
	$(MPISMP)

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libaztecoo_a_AR = $(ALTERNATE_AR)
else

libaztecoo_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libaztecoo_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(AZTECOO_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(AZTECOO_H)

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: AztecOO_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = AztecOO_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






