# @HEADER
# ************************************************************************
# 
#        Claps: A Collection of Domain Decomposition Preconditioners
#               and Solvers
#        Copyright (2006) Sandia Corporation
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
# Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package claps subdirectory src built successfully."
	@echo ""

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
	$(srcdir)/Claps_ConfigDefs.hpp

CORE = \
	$(srcdir)/Claps_Version.cpp

EXTRA_CLAPS_H = \
	$(srcdir)/preconditioners/CLOP_constraint.hpp \
	$(srcdir)/preconditioners/CLOP_graph.hpp \
	$(srcdir)/preconditioners/CLOP_solver.hpp \
	$(srcdir)/preconditioners/CLOP_sub.hpp \
	$(srcdir)/preconditioners/CRD_utils.hpp \
	$(srcdir)/preconditioners/EpetraExtCD_MatrixMatrix.hpp \
	$(srcdir)/preconditioners/my_feti_sparse_solver.hpp \
	$(srcdir)/preconditioners/myzero.hpp \
	$(srcdir)/preconditioners/sparse_lu.hpp \
	$(srcdir)/preconditioners/CLIP_solver.hpp \
	$(srcdir)/preconditioners/sparse_lu2.hpp \
	$(srcdir)/preconditioners/solver_crd.hpp \
	$(srcdir)/preconditioners/CLIP_solver2.hpp \
	$(srcdir)/preconditioners/CRS_serial.hpp \
	$(srcdir)/preconditioners/preconditioner_crd.hpp \
	$(srcdir)/preconditioners/krylov_crd.hpp

EXTRA_CLAPS = \
	$(srcdir)/preconditioners/CLOP_constraint.cpp \
	$(srcdir)/preconditioners/CLOP_graph.cpp \
	$(srcdir)/preconditioners/CLOP_solver.cpp \
	$(srcdir)/preconditioners/CLOP_sub.cpp \
	$(srcdir)/preconditioners/CRD_utils.cpp \
	$(srcdir)/preconditioners/EpetraExtCD_MatrixMatrix.cpp \
	$(srcdir)/preconditioners/sparse_lu.cpp \
	$(srcdir)/preconditioners/CLIP_solver.cpp \
	$(srcdir)/preconditioners/sparse_lu2.cpp \
	$(srcdir)/preconditioners/CLIP_solver2.cpp \
	$(srcdir)/preconditioners/CRS_serial.cpp \
	$(srcdir)/preconditioners/krylov_crd.cpp

EXTRA_SCLAPS_H = \
	$(srcdir)/salinas/SCLOP_solver.hpp \
	$(srcdir)/salinas/SCLIP_solver.hpp \
	$(srcdir)/salinas/CLIP_interface.H \
	$(srcdir)/salinas/CLOP_interface.H \
	$(srcdir)/salinas/MpcLocal.H 

EXTRA_SCLAPS = \
	$(srcdir)/salinas/SCLOP_interface.cpp \
	$(srcdir)/salinas/SCLOP_solver.cpp \
	$(srcdir)/salinas/SCLIP_interface.cpp \
	$(srcdir)/salinas/SCLIP_solver.cpp

if HAVE_SCLAPS

SCLAPS = $(EXTRA_SCLAPS)
SCLAPS_H = $(EXTRA_SCLAPS_H)

else

SCLAPS =
SCLAPS_H =
                                                                                
endif

#
#np#  Files which are included in select builds are typically have
#np#  names that start with EXTRA
#
#np#  EXTRA_MPI is control by --with-mpi and --without-mpi which is
#np#  set up by TAC_ARC_CONFIG_MPI 
#

if HAVE_MPI 
MPI_H = \
	$(EXTRA_CLAPS_H) \
	$(SCLAPS_H)
MPI = \
	$(EXTRA_CLAPS) \
	$(SCLAPS)

else

MPI_H =
MPI =

endif

include $(top_builddir)/Makefile.export.claps

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(CLAPS_INCLUDES))
else
EXPORT_INCLUDES = $(CLAPS_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

#np# The following section can be duplicated for building multiple libraries -
#np# make sure to list the correct categories with each library.
#np# replace new_package with the name of the package being autotool'ed here
# ------------------------------------------------------------------------
# Claps library specifications
# ------------------------------------------------------------------------
CLAPS_LIB = libclaps.a

CLAPS_H = \
	$(CORE_H) \
	$(MPI_H) 

libclaps_a_SOURCES = \
	$(CORE) \
	$(MPI) 

EXTRA_libclaps_a_SOURCES = \
	$(EXTRA_CLAPS) \
	$(EXTRA_CLAPS_H) \
	$(EXTRA_SCLAPS) \
        $(EXTRA_SCLAPS_H) 

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
libclaps_a_AR = $(ALTERNATE_AR)
else

#np# replace new_package with the name of the package being autotool'ed here
libclaps_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
libclaps_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
lib_LIBRARIES = $(CLAPS_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
include_HEADERS = $(CLAPS_H) 

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
nodist_execinclude_HEADERS = Claps_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






