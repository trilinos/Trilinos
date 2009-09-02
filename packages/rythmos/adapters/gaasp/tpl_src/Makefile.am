#@HEADER
# ***********************************************************************
#
#                           Rythmos Package
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
# Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
#
# ***********************************************************************
#@HEADER

# The following line helps the test harness recover from build errors.

all-local:

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

SOURCE_H = \
  $(srcdir)/GAdjointSolve.h\
  $(srcdir)/GErrorEstimate.h\
  $(srcdir)/GForwardSolve.h\
  $(srcdir)/GMeshRefine.h\
  $(srcdir)/GModel.h\
  $(srcdir)/GModelBase.h\
  $(srcdir)/InitGaaspOO.h\
  $(srcdir)/MemoryMgmt.h\
  $(srcdir)/Random.h\
  $(srcdir)/dense.h\
  $(srcdir)/llnlmath.h\
  $(srcdir)/llnltyps.h\
  $(srcdir)/nvector.h

SOURCE = \
  $(srcdir)/AdjointJacobian.cpp \
  $(srcdir)/FDMethod.cpp \
  $(srcdir)/FMethod.cpp \
  $(srcdir)/ForwardJacobian.cpp \
  $(srcdir)/GAdjointSolve.cpp \
  $(srcdir)/GErrorEstimate.cpp \
  $(srcdir)/GForwardSolve.cpp \
  $(srcdir)/GMeshRefine.cpp \
  $(srcdir)/GModel.cpp \
  $(srcdir)/GaaspOO.cpp \
  $(srcdir)/L2Norm.cpp \
  $(srcdir)/MatrixSolve.cpp \
  $(srcdir)/Newton.cpp \
  $(srcdir)/NewtonStep.cpp \
  $(srcdir)/Projections.cpp \
  $(srcdir)/Psi.cpp \
  $(srcdir)/Random.cpp \
  $(srcdir)/Spline.cpp \
  $(srcdir)/dense.cpp \
  $(srcdir)/llnlmath.cpp \
  $(srcdir)/nvector.cpp

#Specify locations of include files within Trilinos

#include $(top_builddir)/Makefile.export.rythmos

#if USING_GNUMAKE
#EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(RYTHMOS_INCLUDES))
#else
#EXPORT_INCLUDES = $(RYTHMOS_INCLUDES)
#endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# ------------------------------------------------------------------------
# GAASP library specifications
# ------------------------------------------------------------------------

GAASP_LIB = libgaasp.a

GAASP_H = \
  $(SOURCE_H) 

libgaasp_a_SOURCES = \
  $(SOURCE) 

EXTRA_libgaasp_a_SOURCES = 

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libgaasp_a_AR = $(ALTERNATE_AR)
else

libgaasp_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libgaasp_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(GAASP_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(GAASP_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include

nodist_execinclude_HEADERS = 

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 
