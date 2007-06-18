# @HEADER
# ************************************************************************
# 
#                 Amesos: Direct Sparse Solver Package
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

SOURCE_H = \
	$(srcdir)/Thyra_AmesosTypes.hpp \
	$(srcdir)/Thyra_AmesosLinearOpWithSolve.hpp \
	$(srcdir)/Thyra_AmesosLinearOpWithSolveFactory.hpp

SOURCE = \
	$(srcdir)/Thyra_AmesosTypes.cpp \
	$(srcdir)/Thyra_AmesosLinearOpWithSolve.cpp \
	$(srcdir)/Thyra_AmesosLinearOpWithSolveFactory.cpp

# Specify locations of include files within Trilinos

include $(top_builddir)/Makefile.export.stratimikos

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(STRATIMIKOS_INCLUDES))
else
EXPORT_INCLUDES = $(STRATIMIKOS_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# --------------------------------------------------------------------
# Amesos Thyra library specifications
# --------------------------------------------------------------------

STRATIMIKOS_AMESOS_LIB = libstratimikosamesos.a

STRATIMIKOS_AMESOS_H = \
	$(top_srcdir)/src/Stratimikos_ConfigDefs.hpp \
	$(SOURCE_H)

libstratimikosamesos_a_SOURCES = \
	$(SOURCE)

#We need the following to support separate build directories

# --------------------------------------------------------------------
# For using a special archiver
# --------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libstratimikosamesos_a_AR = $(ALTERNATE_AR)

else

libstratimikosamesos_a_AR = $(AR) cru

endif

# --------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# --------------------------------------------------------------------

# I have not been able to find any kind of "or" for Automake
# conditionals. JW

 libstratimikosamesos_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# --------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# --------------------------------------------------------------------

lib_LIBRARIES = $(STRATIMIKOS_AMESOS_LIB)

# --------------------------------------------------------------------
# List of all headers to install in $(includedir)
# --------------------------------------------------------------------

include_HEADERS = $(STRATIMIKOS_AMESOS_H) 

# --------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# --------------------------------------------------------------------

execincludedir = $(exec_prefix)/include

# --------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# --------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 
