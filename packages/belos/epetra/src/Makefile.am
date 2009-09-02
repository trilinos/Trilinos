#@HEADER
# ************************************************************************
#
#               Belos: Block Linear Solvers Package
#               Copyright (2001) Sandia Corporation
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
#@HEADER

#np# Note the following paragraph.
# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
#np# There is, at present, no distinction between most of the categories, 
#np# other than those which start with EXTRA_.  All categories which do not
#np# start with EXTRA_ are included in all builds.  Categories which start
#np# with EXTRA_ are included in some builds (typically controlled by 
#np# --with-xyz or --without-xyz)
#
#np# Make sure to list all source files in one of the following categories.
CORE_H = \
	$(srcdir)/BelosEpetraAdapter.hpp \
	$(srcdir)/BelosEpetraOperator.h

CORE = \
	$(srcdir)/BelosEpetraAdapter.cpp \
	$(srcdir)/BelosEpetraOperator.cpp

#
#np#  Files which are included in select builds are typically have
#np#  names that start with EXTRA
#
#np#  EXTRA_MPI is control by --with-mpi and --without-mpi which is
#np#  set up by TAC_ARC_CONFIG_MPI 
#

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

# ------------------------------------------------------------------------
# Includes for required and optional library dependencies
# ------------------------------------------------------------------------
#np# List locations of include files within Trilinos but outside of this 
#np# package
include $(top_builddir)/Makefile.export.belos

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(BELOS_INCLUDES))
else
EXPORT_INCLUDES = $(BELOS_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# ------------------------------------------------------------------------	
#np# The following section can be duplicated for building multiple libraries -
#np# make sure to list the correct categories with each library.
#np# replace new_package with the name of the package being autotool'ed here
# ------------------------------------------------------------------------
# New_Package library specifications
# ------------------------------------------------------------------------
#np# replace new_package with the name of the package being autotool'ed here
BELOS_EPETRA_LIB = libbelosepetra.a

#np# replace new_package with the name of the package being autotool'ed here
BELOS_EPETRA_H = \
	$(top_srcdir)/src/BelosConfigDefs.hpp \
	$(CORE_H) \
	$(MPI_H) 

#np# replace new_package with the name of the package being autotool'ed here
libbelosepetra_a_SOURCES = \
	$(CORE) \
	$(MPI)

#np# replace new_package with the name of the package being autotool'ed here
#EXTRA_libbelosepetra_a_SOURCES = \
#	$(EXTRA_MPI) \
#	$(EXTRA_MPI_H) 

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

#np# replace new_package with the name of the package being autotool'ed here
libbelosepetra_a_AR = $(ALTERNATE_AR)

else

#np# replace new_package with the name of the package being autotool'ed here
libbelosepetra_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
libbelosepetra_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
lib_LIBRARIES = $(BELOS_EPETRA_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
include_HEADERS = $(BELOS_EPETRA_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: New_Package_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 
