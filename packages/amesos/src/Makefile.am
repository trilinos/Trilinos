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


# -----------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside Amesos. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

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


#
#  If the parameter list cannot be made portable to all 
#  platforms, we will have to add a NEED_PARAMETER_LIST
#  flag to configure.ac 
#
if NEED_PARAMETER_LIST
# PARAM_H = \
# 	$(srcdir)/Amesos_Parameter_List.h \
# 	$(srcdir)/Amesos_Parameter_Arbitrary.h \
# 	$(srcdir)/Amesos_Parameter_Entry.h
# PARAM = \
# 	$(srcdir)/Amesos_Parameter_List.cpp \
# 	$(srcdir)/Amesos_Parameter_Entry.cpp
PARAM_H = 
PARAM = 
else
PARAM_H = 
PARAM = 
endif


#
#       We no longer support DscpackOO and hence it is not copied to 
#  the include directory anymore.
#
#	$(srcdir)/DscpackOO.h \
#

if HAVE_AMESOS_DSCPACK
DSCPACK_H = \
	$(srcdir)/Amesos_Dscpack.h
DSCPACK = \
	$(srcdir)/DscpackOO.cpp \
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

if HAVE_AMESOS_SCALAPACK
SCALAPACK_H = \
	$(srcdir)/Amesos_Scalapack.h
SCALAPACK = \
	$(srcdir)/Amesos_Scalapack.cpp
else
SCALAPACK_H =
SCALAPACK =
endif

#
#  Here is where we will add Tim Davis' files 
#
if HAVE_AMESOS_KLU
#
#  These files do not need to be copied to the include directory on install,
#  hence they are commented out from here.  
#
# KLU_H = \
# 	$(srcdir)/klu/klu_btf.h \
# 	$(srcdir)/Amesos_Klu.h \
# 	$(srcdir)/klu/cbtf.h \
# 	$(srcdir)/klu/klu.h \
# 	$(srcdir)/klu/klu_dump.h \
# 	$(srcdir)/klu/klu_kernel.h
KLU_H = \
	$(srcdir)/Amesos_Klu.h 

KLU = \
	$(srcdir)/Amesos_Klu.cpp \
	$(srcdir)/AMD/Source/amd_1.c \
	$(srcdir)/AMD/Source/amd_2.c \
	$(srcdir)/AMD/Source/amd_aat.c \
	$(srcdir)/AMD/Source/amd_control.c \
	$(srcdir)/AMD/Source/amd_defaults.c \
	$(srcdir)/AMD/Source/amd_dump.c \
	$(srcdir)/AMD/Source/amd_info.c \
	$(srcdir)/AMD/Source/amd_order.c \
	$(srcdir)/AMD/Source/amd_postorder.c \
	$(srcdir)/AMD/Source/amd_post_tree.c \
	$(srcdir)/AMD/Source/amd_preprocess.c \
	$(srcdir)/AMD/Source/amd_valid.c \
	$(srcdir)/btf/btf.f \
	$(srcdir)/btf/cbtf.c \
	$(srcdir)/btf/mc13d.f \
	$(srcdir)/btf/mc13e.f \
	$(srcdir)/btf/mc21a.f \
	$(srcdir)/btf/mc21b.f \
	$(srcdir)/colamd/colamd.c \
	$(srcdir)/klu_btf/klu_btf_analyze.c \
	$(srcdir)/klu_btf/klu_btf_defaults.c \
	$(srcdir)/klu_btf/klu_btf_factor.c \
	$(srcdir)/klu_btf/klu_btf_free_numeric.c \
	$(srcdir)/klu_btf/klu_btf_free_symbolic.c \
	$(srcdir)/klu_btf/klu_btf_refactor.c \
	$(srcdir)/klu_btf/klu_btf_scale.c \
	$(srcdir)/klu_btf/klu_btf_solve.c \
	$(srcdir)/klu_btf/klu_dump.c \
	$(srcdir)/klu_kernel/klu.c \
	$(srcdir)/klu_kernel/klu_kernel.c

else
KLU_H =
KLU =
endif

if HAVE_AMESOS_MUMPS
MUMPS_H = \
	$(srcdir)/Amesos_Mumps.h \
	$(srcdir)/Amesos_EpetraBaseSolver.h \
	$(srcdir)/EpetraExt_Redistor.h
MUMPS = \
	$(srcdir)/Amesos_Mumps.cpp \
	$(srcdir)/Amesos_EpetraBaseSolver.cpp \
	$(srcdir)/EpetraExt_Redistor.cpp
else
MUMPS_H =
MUMPS =
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
       $(srcdir)/ksparse/driver.c \
       $(srcdir)/ksparse/spalloc.c \
       $(srcdir)/ksparse/spfactor.c \
       $(srcdir)/ksparse/spsmp.c \
       $(srcdir)/ksparse/sputils.c
else
KUNDERT_H = 
KUNDERT =
endif

#
#  This file does not need to be copied to the include directory
#	$(srcdir)/CrsMatrixTranspose.h
COMP_H = \
	Amesos_BaseSolver.h \
	Amesos_ConfigDefs.h 
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

AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src \
	-I$(top_srcdir)/src/AMD/Include \
	-I$(top_srcdir)/src/AMD/Source \
	-I$(top_srcdir)/src/colamd \
	-I$(top_srcdir)/src/btf \
	-I$(top_srcdir)/src/klu_btf \
	-I$(top_srcdir)/src/klu_kernel \
		-I$(top_srcdir)/../teuchos/src  -I$(top_builddir)/../teuchos/src \
		-I$(top_builddir)/../epetra/src \
	      $(SLUSERIAL_CPPFLAGS) $(KUNDERT_CPPFLAGS)

# ------------------------------------------------------------------------
# Amesos library specifications
# ------------------------------------------------------------------------

AMESOS_LIB = libamesos.a


AMESOS_H = \
	$(PARAM_H) \
	$(DSCPACK_H) \
	$(UMFPACK_H) \
	$(SCALAPACK_H) \
	$(KLU_H) \
	$(MUMPS_H) \
	$(FACTORY_H) \
	$(SUPERLU_H) \
	$(SUPERLUDIST_H) \
	$(KUNDERT_H) \
	$(COMP_H) \
	$(srcdir)/Amesos_BaseSolver.h


libamesos_a_SOURCES = \
	$(PARAM) \
	$(DSCPACK) \
	$(KLU) \
	$(UMFPACK) \
	$(SCALAPACK) \
	$(MUMPS) \
	$(FACTORY) \
	$(SUPERLU) \
	$(SUPERLUDIST) \
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






