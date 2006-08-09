# @HEADER
# ************************************************************************
#
#            Galeri: Finite Element and Matrix Generation Package
#                 Copyright (2006) ETHZ/Sandia Corporation
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
#
# Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
#
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package galeri subdirectory src-pfem built successfully."
	@echo ""

CORE_H = $(srcdir)/Galeri_problem_ScalarLaplacian.h

CORE = $(srcdir)/Galeri_core_Utils.cpp \
       $(srcdir)/Galeri_core_Constants.cpp \
       $(srcdir)/Galeri_grid_Generator.cpp

include $(top_builddir)/Makefile.export.galeri

AM_CPPFLAGS = $(GALERI_INCLUDES)

# ------------------------------------------------------------------------
# pfem library specifications
# ------------------------------------------------------------------------
PFEM_LIB = libpfem.a

PFEM_H = $(CORE_H) 

libpfem_a_SOURCES = $(CORE) 

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR
libpfem_a_AR = $(ALTERNATE_AR)
else
libpfem_a_AR = $(AR) cru
endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libpfem_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(PFEM_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(PFEM_H) 

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
	$(MAKE) libpfem.a
