# @HEADER
# ************************************************************************
# 
#               Epetra: Linear Algebra Services Package 
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
# might someday be needed by some other code outside Epetra. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

CORE = \
	$(srcdir)/Epetra_test_functions.cpp \
	$(srcdir)/Epetra_matrix_data.cpp

#Although these headers do not have to be installed, they do need to be
#included in the distribution tarball.
EXTRA_DIST = \
	$(srcdir)/Epetra_test_functions.h \
	$(srcdir)/Epetra_matrix_data.h

DEFAULT_INCLUDES = \
	-I$(top_srcdir)/test/src \
	-I$(top_builddir)/src \
	-I$(top_srcdir)/src

# ------------------------------------------------------------------------
# Epetra_test library specifications
# ------------------------------------------------------------------------

EPETRA_TEST_LIB = libepetra_test.a

libepetra_test_a_SOURCES = \
	$(CORE)


#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libepetra_test_a_AR = $(ALTERNATE_AR)
else

libepetra_test_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# On alpha-dec machines, the extra .o's are needed for Epetra only.
# ------------------------------------------------------------------------

libepetra_test_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

lib_LIBRARIES = $(EPETRA_TEST_LIB)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

