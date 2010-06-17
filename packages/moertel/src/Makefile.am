# @HEADER
# ************************************************************************
# 
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
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
# Questions? Contact Glen Hansen (Glen.Hansen@inl.gov)
# 
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package moertel subdirectory src built successfully."
	@echo ""

CORE = \
	$(srcdir)/mrtr_convexhull.cpp \
	$(srcdir)/mrtr_point.cpp \
	$(srcdir)/mrtr_functions.cpp \
	$(srcdir)/mrtr_projector3D.cpp \
	$(srcdir)/mrtr_integrator.cpp \
	$(srcdir)/mrtr_projector.cpp \
	$(srcdir)/mrtr_interface_complete.cpp \
	$(srcdir)/mrtr_segment_bilinearquad.cpp \
	$(srcdir)/mrtr_interface_integrate3D.cpp \
	$(srcdir)/mrtr_segment_bilineartri.cpp \
	$(srcdir)/mrtr_interface_integrate.cpp \
	$(srcdir)/mrtr_interface_project.cpp \
	$(srcdir)/mrtr_segment.cpp \
	$(srcdir)/mrtr_interface_tools.cpp \
	$(srcdir)/mrtr_manager.cpp \
	$(srcdir)/mrtr_segment_linear1D.cpp \
	$(srcdir)/mrtr_manager_solve.cpp \
	$(srcdir)/mrtr_solver.cpp \
	$(srcdir)/mrtr_ml_preconditioner.cpp \
	$(srcdir)/mrtr_node.cpp \
	$(srcdir)/mrtr_overlap.cpp \
	$(srcdir)/mrtr_overlap_utils.cpp \
	$(srcdir)/mrtr_pnode.cpp


CORE_H = \
	$(srcdir)/mrtr_pnode.H \
	$(srcdir)/mrtr_function.H \
	$(srcdir)/mrtr_point.H \
	$(srcdir)/mrtr_functions.H \
	$(srcdir)/mrtr_integrator.H \
	$(srcdir)/mrtr_projector.H \
	$(srcdir)/mrtr_interface.H \
	$(srcdir)/mrtr_segment_bilinearquad.H \
	$(srcdir)/mrtr_segment_bilineartri.H \
	$(srcdir)/mrtr_segment.H \
	$(srcdir)/mrtr_manager.H \
	$(srcdir)/mrtr_segment_linear1D.H \
	$(srcdir)/mrtr_solver.H \
	$(srcdir)/mrtr_ml_preconditioner.H \
	$(srcdir)/mrtr_node.H \
	$(srcdir)/mrtr_overlap.H

UTIL_H = \
	$(srcdir)/mrtr_utils.H

UTIL = \
	$(srcdir)/mrtr_utils.cpp 

include $(top_builddir)/Makefile.export.moertel

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(MOERTEL_INCLUDES))
else
EXPORT_INCLUDES = $(MOERTEL_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# ------------------------------------------------------------------------
# Moertel library specifications
# ------------------------------------------------------------------------

MOERTEL_LIB = libmoertel.a

MOERTEL_H = \
	$(CORE_H) \
	$(UTIL_H)

libmoertel_a_SOURCES = \
	$(CORE) \
	$(UTIL)

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libmoertel_a_AR = $(ALTERNATE_AR)
else

libmoertel_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libmoertel_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(MOERTEL_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(MOERTEL_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Moertel_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(includedir)
nodist_execinclude_HEADERS = Moertel_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

