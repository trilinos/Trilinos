# @HEADER
# ************************************************************************
# 
#          Tpetra: Templated Linear Algebra Services Package
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

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

CORE_H = \
	$(srcdir)/Tpetra_Object.hpp \
	$(srcdir)/Tpetra_Platform.hpp \
	$(srcdir)/Tpetra_Comm.hpp \
	$(srcdir)/Tpetra_Distributor.hpp \
	$(srcdir)/Tpetra_Directory.hpp

CORE = 

UTIL_H = \
	$(srcdir)/Tpetra_ConfigDefs.hpp \
	$(srcdir)/Tpetra_PacketTraits.hpp \
	$(srcdir)/Tpetra_Version.hpp \
	$(srcdir)/Tpetra_Util.hpp

UTIL = \
	$(srcdir)/Tpetra_Version.cpp

ENUM_H = \
	$(srcdir)/Tpetra_CombineMode.hpp

UNIPROC_H = 

UNIPROC = 

SERIAL_H = \
	$(srcdir)/Tpetra_SerialPlatform.hpp \
	$(srcdir)/Tpetra_SerialComm.hpp \
	$(srcdir)/Tpetra_SerialDistributor.hpp 

SERIAL = 

EXTRA_MPI_H = \
	$(srcdir)/Tpetra_MpiPlatform.hpp \
	$(srcdir)/Tpetra_MpiComm.hpp \
	$(srcdir)/Tpetra_MpiDistributor.hpp \
	$(srcdir)/Tpetra_MpiData.hpp \
	$(srcdir)/Tpetra_MpiTraits.hpp 

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

DIST_H = \
	$(srcdir)/Tpetra_ElementSpace.hpp \
	$(srcdir)/Tpetra_ElementSpaceData.hpp \
	$(srcdir)/Tpetra_BlockElementSpace.hpp \
	$(srcdir)/Tpetra_BlockElementSpaceData.hpp \
	$(srcdir)/Tpetra_VectorSpace.hpp \
	$(srcdir)/Tpetra_VectorSpaceData.hpp \
	$(srcdir)/Tpetra_Import.hpp \
	$(srcdir)/Tpetra_ImportData.hpp \
	$(srcdir)/Tpetra_BasicDirectory.hpp

DIST = 

COMP_H = \
	$(srcdir)/Tpetra_Vector.hpp \
	$(srcdir)/Tpetra_VectorData.hpp \
	$(srcdir)/Tpetra_CisMatrix.hpp \
	$(srcdir)/Tpetra_CisMatrixData.hpp

COMP = 


# List locations of include files within Trilinos but outside of this package

AM_CPPFLAGS = -I$(top_srcdir)/../teuchos/src -I$(top_builddir)/../teuchos/src

# The following section can be duplicated for building multiple libraries -
# make sure to list the correct categories with each library.
# ------------------------------------------------------------------------
# Tpetra library specifications
# ------------------------------------------------------------------------
TPETRA_LIB = libtpetra.a

TPETRA_H = \
	$(CORE_H) \
	$(UTIL_H) \
	$(ENUM_H) \
	$(UNIPROC_H) \
	$(SERIAL_H) \
	$(USERBASE_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H) 

libtpetra_a_SOURCES = \
	$(CORE) \
	$(UTIL) \
	$(UNIPROC) \
	$(SERIAL) \
	$(DIST) \
	$(COMP) \
	$(MPI)

EXTRA_libtpetra_a_SOURCES = \
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

libtpetra_a_AR = $(ALTERNATE_AR)
else

libtpetra_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libtpetra_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(TPETRA_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(TPETRA_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Tpetra_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = Tpetra_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






