# @HEADER
# ************************************************************************
# 
#                 TriUtils: Trilinos Utilities Package
#                 Copyright (2001) Sandia Corporation
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
# For now, these are just guesses as to where the files belong.
# Files that prove to be troublesome have been commented out and 
# placed at the bottom of the appropriate section

CORE_H = \
	$(srcdir)/Trilinos_Util.h \
	$(srcdir)/Trilinos_Util_ReadMatrixMarket2Epetra.h \
	$(srcdir)/Trilinos_Util_ReadTriples2Epetra.h \
	$(srcdir)/Trilinos_Util_CountMatrixMarket.h \
	$(srcdir)/Trilinos_Util_CountTriples.h \
	$(srcdir)/Trilinos_Util_CrsMatrixGallery.h \
	$(srcdir)/Trilinos_Util_CommandLineParser.h

CORE = \
	$(srcdir)/Trilinos_Util_duscr_vbr.cpp \
	$(srcdir)/Trilinos_Util_dusds_vbr.cpp \
	$(srcdir)/Trilinos_Util_distrib_vbr_matrix.cpp \
	$(srcdir)/Trilinos_Util_read_hb.cpp \
	$(srcdir)/Trilinos_Util_svbrres.cpp \
	$(srcdir)/Trilinos_Util_dusmm.cpp \
	$(srcdir)/Trilinos_Util_scscmv.cpp \
	$(srcdir)/Trilinos_Util_read_vec.cpp \
	$(srcdir)/Trilinos_Util_write_vec.cpp \
	$(srcdir)/Trilinos_Util_create_vbr.cpp \
	$(srcdir)/Trilinos_Util_msr2vbr.cpp \
	$(srcdir)/Trilinos_Util_scscres.cpp \
	$(srcdir)/Trilinos_Util_distrib_msr_matrix.cpp \
	$(srcdir)/Trilinos_Util_read_coo.cpp \
	$(srcdir)/Trilinos_Util_smsrres.cpp \
	$(srcdir)/Trilinos_Util_csrcsc.cpp \
	$(srcdir)/Trilinos_Util_csrmsr.cpp \
	$(srcdir)/Trilinos_Util_ssrcsr.cpp \
	$(srcdir)/Trilinos_Util_coocsr.cpp \
	$(srcdir)/Trilinos_Util_ReadMatrixMarket2Epetra.cpp \
	$(srcdir)/Trilinos_Util_ReadTriples2Epetra.cpp \
	$(srcdir)/Trilinos_Util_CountMatrixMarket.cpp \
	$(srcdir)/Trilinos_Util_CountTriples.cpp \
	$(srcdir)/Trilinos_Util_ReadHb2Epetra.cpp \
	$(srcdir)/Trilinos_Util_ReadHpc2Epetra.cpp \
	$(srcdir)/Trilinos_Util_ReadHb2EpetraVbr.cpp \
	$(srcdir)/Trilinos_Util_GenerateCrsProblem.cpp \
	$(srcdir)/Trilinos_Util_GenerateVbrProblem.cpp \
	$(srcdir)/Trilinos_Util_CommandLineParser.cpp \
	$(srcdir)/Trilinos_Util_CrsMatrixGallery.cpp


AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src

WITH_CPP_H =

WITH_CPP =


UTIL_H = \
	$(srcdir)/iohb.h \
	$(srcdir)/Trilinos_Util_Version.h

UTIL = \
	$(srcdir)/iohb.cpp

ENUM_H = 

UNIPROC_H = 

UNIPROC = 

MPI_H = 
MPI = 

SERIAL_H = 
SERIAL = 

#Triutils does not support threads
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

LANG_WRAPS = 

DIST_H = 

DIST = 

COMP_H = 

COMP = 

# ------------------------------------------------------------------------
# Aztec library specifications
# ------------------------------------------------------------------------

TRIUTILS_LIB = libtriutils.a


TRIUTILS_H = \
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

libtriutils_a_SOURCES = \
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

libtriutils_a_AR = $(ALTERNATE_AR)
else

libtriutils_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libtriutils_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(TRIUTILS_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(TRIUTILS_H)

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Triutils_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = Triutils_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






