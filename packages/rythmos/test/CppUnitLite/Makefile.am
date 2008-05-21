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

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

CORE_H = \
  $(srcdir)/TestHarness.h\
  $(srcdir)/Test.h\
  $(srcdir)/TestRegistry.h\
  $(srcdir)/TestResult.h\
  $(srcdir)/SimpleString.h\
  $(srcdir)/Failure.h 

CORE = \
  $(srcdir)/Test.cpp\
  $(srcdir)/TestRegistry.cpp\
  $(srcdir)/TestResult.cpp\
  $(srcdir)/SimpleString.cpp\
  $(srcdir)/Failure.cpp

UTIL_H =

UTIL = 

ENUM_H = 

UNIPROC_H = 

UNIPROC = 

SERIAL_H = 

SERIAL = 

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

USERBASE_H = 


LANG_WRAPS_H = 

LANG_WRAPS = 

DIST_H = 

DIST = 

COMP_H =

COMP = 

#Specify locations of include files within Trilinos

include $(top_builddir)/Makefile.export.rythmos

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(RYTHMOS_INCLUDES))
else
EXPORT_INCLUDES = $(RYTHMOS_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# ------------------------------------------------------------------------
# Rythmos library specifications
# ------------------------------------------------------------------------

CPPUNITLITE_LIB = libcppunitlite.a

CPPUNITLITE_H = \
  $(CORE_H) \
  $(UTIL_H) \
  $(ENUM_H) \
  $(UNIPROC_H) \
  $(SERIAL_H) \
  $(USERBASE_H) \
  $(DIST_H) \
  $(COMP_H) \
  $(MPI_H)

libcppunitlite_a_SOURCES = \
  $(CORE) \
  $(UTIL) \
  $(UNIPROC) \
  $(SERIAL) \
  $(DIST) \
  $(COMP) \
  $(MPI)

EXTRA_libcppunitlite_a_SOURCES = \
  $(EXTRA_MPI) \
  $(EXTRA_MPI_H)

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libcppunitlite_a_AR = $(ALTERNATE_AR)
else

libcppunitlite_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libcppunitlite_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(CPPUNITLITE_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(CPPUNITLITE_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include

#nodist_execinclude_HEADERS = Rythmos_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 
