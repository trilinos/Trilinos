# @HEADER
# ************************************************************************
# 
#     EpetraExt: Epetra Extended - Linear Algebra Services Package 
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
# might someday be needed by some other code outside EpetraExt. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

INOUT_H = \
	$(srcdir)/EpetraExt_mmio.h \
     $(srcdir)/EpetraExt_RowMatrixOut.h \
     $(srcdir)/EpetraExt_MultiVectorOut.h \
     $(srcdir)/EpetraExt_VectorOut.h \
     $(srcdir)/EpetraExt_BlockMapOut.h \
     $(srcdir)/EpetraExt_BlockMapIn.h \
     $(srcdir)/EpetraExt_CrsMatrixIn.h \
     $(srcdir)/EpetraExt_MultiVectorIn.h \
     $(srcdir)/EpetraExt_VectorIn.h \
     $(srcdir)/EpetraExt_Exception.h \
     $(srcdir)/EpetraExt_Utils.h \
     $(srcdir)/EpetraExt_HDF5.h

INOUT = \
	$(srcdir)/EpetraExt_mmio.cpp \
     $(srcdir)/EpetraExt_RowMatrixOut.cpp \
     $(srcdir)/EpetraExt_MultiVectorOut.cpp \
     $(srcdir)/EpetraExt_VectorOut.cpp \
     $(srcdir)/EpetraExt_BlockMapOut.cpp \
     $(srcdir)/EpetraExt_BlockMapIn.cpp \
     $(srcdir)/EpetraExt_CrsMatrixIn.cpp \
     $(srcdir)/EpetraExt_MultiVectorIn.cpp \
     $(srcdir)/EpetraExt_VectorIn.cpp \
     $(srcdir)/EpetraExt_Exception.cpp \
     $(srcdir)/EpetraExt_Utils.cpp \
     $(srcdir)/EpetraExt_HDF5.cpp

# ------------------------------------------------------------------------
# EpetraExt inout library specifications
# ------------------------------------------------------------------------

libepetraext_inout_a_SOURCES = \
	$(INOUT)

#If these are not listed, the files will not be included in
#the distribution tarball.
EXTRA_libepetraext_inout_a_SOURCES = 

if BUILD_INOUT

# libraries to be installed in $(libexecdir)
EPETRAEXT_INOUT_LIB = libepetraext_inout.a
EPETRAEXT_INOUT_H = \
	$(INOUT_H)

else

EPETRAEXT_INOUT_LIB =
EPETRAEXT_INOUT_H =

endif


#We need the following to support separate build directories
AM_CPPFLAGS=-I$(top_builddir)/../epetra/src -I$(top_srcdir)/../epetra/src

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libepetraext_inout_a_AR = $(ALTERNATE_AR)

else

libepetraext_inout_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# On alpha-dec machines, the extra .o's are needed for Epetra only.
# ------------------------------------------------------------------------

libepetraext_inout_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(EPETRAEXT_INOUT_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(EPETRAEXT_INOUT_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Epetra_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS =../EpetraExt_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Files to be deleted by 'make clean'
# ------------------------------------------------------------------------

MOSTLYCLEANFILES = $(TRILINOS_TEMPLATE_OBJS)

