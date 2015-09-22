# @HEADER
# ************************************************************************
# 
#               Thyra: Trilinos Solver Framework Core
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
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
# There is, at present, no distinction between most of the categories, 
# other than those which start with EXTRA_.  All categories which do not
# start with EXTRA_ are included in all builds.  Categories which start
# with EXTRA_ are included in some builds (typically controlled by 
# --with-xyz or --without-xyz)
#


SOURCE_H = \
  $(srcdir)/Rythmos_GAASPErrorEstimate.hpp \
  $(srcdir)/Rythmos_GAASPErrorEstimator.hpp \
  $(srcdir)/Rythmos_GAASPGModel_ThyraModelEvaluator.hpp \
  $(srcdir)/Rythmos_GAASPHelpers.hpp \
  $(srcdir)/Rythmos_GAASPInterface.hpp

SOURCE = \
  $(srcdir)/Rythmos_GAASPErrorEstimate.cpp \
  $(srcdir)/Rythmos_GAASPErrorEstimator.cpp \
  $(srcdir)/Rythmos_GAASPGModel_ThyraModelEvaluator.cpp \
  $(srcdir)/Rythmos_GAASPHelpers.cpp \
  $(srcdir)/Rythmos_GAASPInterface.cpp


# Specify locations of include files within Trilinos

include $(top_builddir)/Makefile.export.rythmos

_EXPORT_INCLUDES = \
  $(GAASP_INCLUDES) \
  $(RYTHMOS_INCLUDES) 

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(_EXPORT_INCLUDES))
else
EXPORT_INCLUDES = $(_EXPORT_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# --------------------------------------------------------------------
# GAASP Rythmos library specifications
# --------------------------------------------------------------------

RYTHMOS_GAASP_LIB = librythmosgaasp.a

EPETRA_THYRA_H = \
	$(top_srcdir)/src/Rythmos_ConfigDefs.hpp \
	$(SOURCE_H)

librythmosgaasp_a_SOURCES = \
	$(SOURCE)

#We need the following to support separate build directories

# --------------------------------------------------------------------
# For using a special archiver
# --------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

librythmosgaasp_a_AR = $(ALTERNATE_AR)

else

librythmosgaasp_a_AR = $(AR) cru

endif

# --------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# --------------------------------------------------------------------

# I have not been able to find any kind of "or" for Automake
# conditionals. JW

librythmosgaasp_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# --------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# --------------------------------------------------------------------

lib_LIBRARIES = $(GAASP_LIBS) $(RYTHMOS_GAASP_LIB)

# --------------------------------------------------------------------
# List of all headers to install in $(includedir)
# --------------------------------------------------------------------

include_HEADERS = $(RYTHMOS_GAASP_H) 

# --------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# --------------------------------------------------------------------

execincludedir = $(exec_prefix)/include

# --------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# --------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 
