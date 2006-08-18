# @HEADER
# ************************************************************************
# 
#             PyTrilinos Python-Trilinos Utility Package
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
	@echo "Trilinos package pytrilinos subdirectory src built successfully."
	@echo ""

CORE_H = \
	$(srcdir)/Callback.h        \
	$(srcdir)/FileStream.h      \
	$(srcdir)/numpy_include.h \
	$(srcdir)/NumPyWrapper.h    \
	$(srcdir)/NumPyArray.h

CORE = \
	$(srcdir)/Callback.cpp     \
	$(srcdir)/FileStream.cpp   \
	$(srcdir)/NumPyWrapper.cpp \
	$(srcdir)/NumPyArray.cpp

EXTRA_DIST = PyEpetra_Utils.h  \
	     PyTeuchos_Utils.h \
	     numpy.i

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


# Include files outside of this package
NUMPY_INC   = $(shell $(PYTHON) -c "from numpy import get_numpy_include; print '-I' + get_numpy_include()")
AM_CPPFLAGS = $(PYTHON_CSPEC) $(NUMPY_INC)

# ------------------------------------------------------------------------
# PyTrilinos library specifications
# ------------------------------------------------------------------------
PYTRILINOS_LIB = libpytrilinos.a

PYTRILINOS_H = \
	$(CORE_H) \
	$(MPI_H)

libpytrilinos_a_SOURCES = \
	$(CORE) \
	$(MPI)

EXTRA_libpytrilinos_a_SOURCES = \
	$(EXTRA_MPI)

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.
if USE_ALTERNATE_AR
  libpytrilinos_a_AR = $(ALTERNATE_AR)
else
  libpytrilinos_a_AR = $(AR) cru
endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------
libpytrilinos_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------
lib_LIBRARIES = $(PYTRILINOS_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------
include_HEADERS = $(PYTRILINOS_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------
# SPECIAL NOTE: Pytrilinos_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = PyTrilinos_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------
MAINTAINERCLEANFILES = Makefile.in 
