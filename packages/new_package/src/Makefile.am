# @HEADER
# ************************************************************************
# 
#                     New_Package Example Package
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
	$(srcdir)/Newp_Hello.h \
	$(srcdir)/New_Package_ConfigDefs.h

CORE = \
	$(srcdir)/Newp_Hello.cpp

UTIL_H = \
	$(srcdir)/New_Package_Version.h

UTIL = 

ENUM_H = 

UNIPROC_H = 

UNIPROC = 

SERIAL_H = 

SERIAL = 

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
#np# This section should come out, or be modified to reflect options in the 
#np# package being autotool'ed
EXTRA_NEWP_SWAHILI_H = \
	$(srcdir)/Newp_Jambo.h 

EXTRA_NEWP_SWAHILI = \
	$(srcdir)/Newp_Jambo.cpp

if HAVE_NEWP_SWAHILI

NEWP_SWAHILI_H = \
	$(EXTRA_NEWP_SWAHILI_H)
NEWP_SWAHILI = \
	$(EXTRA_NEWP_SWAHILI)

else

NEWP_SWAHILI_H =
NEWP_SWAHILI =

endif



USERBASE_H = 


LANG_WRAPS_H = 

LANG_WRAPS = 

DIST_H = 

DIST = 

COMP_H =

COMP = 


#np# List locations of include files within Trilinos but outside of this 
#np# package
#np# - Add packages that this package is dependent on
#np#   Remove epetra if this package does not need it
#  
AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src

#np# The following section can be duplicated for building multiple libraries -
#np# make sure to list the correct categories with each library.
#np# replace new_package with the name of the package being autotool'ed here
# ------------------------------------------------------------------------
# New_Package library specifications
# ------------------------------------------------------------------------
#np# replace new_package with the name of the package being autotool'ed here
NEW_PACKAGE_LIB = libnew_package.a

#np# replace new_package with the name of the package being autotool'ed here
NEW_PACKAGE_H = \
	$(CORE_H) \
	$(UTIL_H) \
	$(ENUM_H) \
	$(UNIPROC_H) \
	$(SERIAL_H) \
	$(USERBASE_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H) \
	$(NEWP_SWAHILI_H)

#np# replace new_package with the name of the package being autotool'ed here
libnew_package_a_SOURCES = \
	$(CORE) \
	$(UTIL) \
	$(UNIPROC) \
	$(SERIAL) \
	$(DIST) \
	$(COMP) \
	$(MPI) \
	$(NEWP_SWAHILI)

#np# replace new_package with the name of the package being autotool'ed here
EXTRA_libnew_package_a_SOURCES = \
	$(EXTRA_MPI) \
	$(EXTRA_MPI_H) \
	$(EXTRA_NEWP_SWAHILI) \
	$(EXTRA_NEWP_SWAHILI_H)

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
libnew_package_a_AR = $(ALTERNATE_AR)
else

#np# replace new_package with the name of the package being autotool'ed here
libnew_package_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
libnew_package_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
lib_LIBRARIES = $(NEW_PACKAGE_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
include_HEADERS = $(NEW_PACKAGE_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
# SPECIAL NOTE: New_Package_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include
#np# replace new_package with the name of the package being autotool'ed here
nodist_execinclude_HEADERS = New_Package_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






