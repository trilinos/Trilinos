# @HEADER
# ************************************************************************
# 
#                     Galeri Example Package
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
# Questions? Contact Marzio Sala (marzio.sala@gmail.com) 
# 
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package galeri subdirectory src-phoenix built successfully."
	@echo ""

CORE_H = $(srcdir)/phx_problem_ScalarLaplacian.h

CORE = $(srcdir)/phx_core_Utils.cpp \
       $(srcdir)/phx_core_Constants.cpp \
       $(srcdir)/phx_grid_Generator.cpp

AM_CPPFLAGS = -I$(srcdir)                      \
              -I$(top_srcdir)/../epetra/src    \
	      -I$(top_builddir)/../epetra/src  \
	      -I$(top_srcdir)/../teuchos/src   \
	      -I$(top_builddir)/../teuchos/src \
	      -I$(top_srcdir)/../epetraext/src/inout   \
	      -I$(top_srcdir)/../epetraext/src   \
	      -I$(top_builddir)/../epetraext/src \
              -I$(srcdir)/../src                 \
              -I$(top_srcdir)/../src    

# ------------------------------------------------------------------------
# phoenix library specifications
# ------------------------------------------------------------------------
PHOENIX_LIB = libphoenix.a

PHOENIX_H = $(CORE_H) 

libphoenix_a_SOURCES = $(CORE) 

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR
libphoenix_a_AR = $(ALTERNATE_AR)
else
libphoenix_a_AR = $(AR) cru
endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libphoenix_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(PHOENIX_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(PHOENIX_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Library Target for Circular Dependencies
# ------------------------------------------------------------------------

all-libs:
	cd $(top_builddir)/../epetra ; $(MAKE) all-libs
	cd $(top_builddir)/../teuchos ; $(MAKE) all-libs
	$(MAKE) libphoenix.a
