# @HEADER
#  ************************************************************************
#  
#                Kokkos: Fast Core Kernels Package
#               Copyright (2004) Sandia Corporation
#  
#  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
#  license for use of this work by or on behalf of the U.S. Government.
#  
#  This library is free software; you can redistribute it and or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation; either version 2.1 of the
#  License, or (at your option) any later version.
#   
#  This library is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#   
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
#  USA
#  Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
#  
#  ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

CORE_H = \
	$(srcdir)/Kokkos_ConfigDefs.hpp \
	$(srcdir)/Kokkos_BaseSparseMultiply.hpp \
     $(srcdir)/Kokkos_BaseSparseSolve.hpp \
	$(srcdir)/Kokkos_CisMatrix.hpp \
	$(srcdir)/Kokkos_DenseMultiVector.hpp \
	$(srcdir)/Kokkos_DenseVector.hpp \
	$(srcdir)/Kokkos_HbMatrix.hpp \
	$(srcdir)/Kokkos_MultiVector.hpp \
	$(srcdir)/Kokkos_PackedSparseMultiply.hpp \
	$(srcdir)/Kokkos_Permutation.hpp \
	$(srcdir)/Kokkos_SparseOperation.hpp \
	$(srcdir)/Kokkos_Vector.hpp

CORE = 

UTIL_H = \
	$(srcdir)/Kokkos_Flops.hpp \
	$(srcdir)/Kokkos_Time.hpp \
	$(srcdir)/Kokkos_Fortran_wrappers.hpp \
	$(srcdir)/Kokkos_f77func.hpp 

UTIL = \
	$(srcdir)/Kokkos_Version.cpp \
	$(srcdir)/Kokkos_dcrsmv.f \
	$(srcdir)/Kokkos_dcrsmm.f 

ENUM_H = 

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

LANG_WRAPS_H = 

LANG_WRAPS = 

DIST_H = 

DIST = 

COMP_H = \
	$(srcdir)/Kokkos_CompObject.hpp

COMP = 


AM_CPPFLAGS =

# ------------------------------------------------------------------------
# Kokkos library specifications
# ------------------------------------------------------------------------

KOKKOS_LIB = libkokkos.a

KOKKOS_H = \
	$(CORE_H) \
	$(UTIL_H) \
	$(ENUM_H) \
	$(SERIAL_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H)

libkokkos_a_SOURCES = \
	$(CORE) \
	$(UTIL) \
	$(SERIAL) \
	$(DIST) \
	$(COMP) \
	$(MPI)

EXTRA_libkokkos_a_SOURCES = \
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

#np# replace new_package with the name of the package being autotool'ed here
libkokkos_a_AR = $(ALTERNATE_AR)
else

#np# replace new_package with the name of the package being autotool'ed here
libkokkos_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
libkokkos_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
lib_LIBRARIES = $(KOKKOS_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
include_HEADERS = $(KOKKOS_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Kokkos_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include
#np# replace new_package with the name of the package being autotool'ed here
nodist_execinclude_HEADERS = Kokkos_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

