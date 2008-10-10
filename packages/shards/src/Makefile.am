# @HEADER
# ************************************************************************
# 
#                 Shards : Shared Discretization Tools
#                  Copyright (2008) Sandia Corporation
# 
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
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
# 
# Questions? Contact Carter Edwards (hcedwar@sandia.gov), Sandia
# National Laboratories.
# 
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package shards subdirectory src built successfully."
	@echo ""


# ------------------------------------------------------------------------
# List of all files
# ------------------------------------------------------------------------

SHARDS_H = \
	$(srcdir)/Shards_CellTopologyData.h

SHARDS_HPP = \
	$(srcdir)/Shards_Array.hpp	\
	$(srcdir)/Shards_CellTopologyTraits.hpp	\
	$(srcdir)/Shards_BasicTopologies.hpp	\
	$(srcdir)/Shards_IndexList.hpp	\
	$(srcdir)/Shards_SimpleArrayOps.hpp	\
	$(srcdir)/Shards_CellTopology.hpp	\
	$(srcdir)/Shards_TypeList.hpp

SHARDS_CPP = \
	$(srcdir)/Shards_Array.cpp \
	$(srcdir)/Shards_BasicTopologies.cpp

# ------------------------------------------------------------------------
# Includes for required and optional library dependencies
# ------------------------------------------------------------------------
include $(top_builddir)/Makefile.export.shards

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(SHARDS_INCLUDES))
else
EXPORT_INCLUDES = $(SHARDS_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

#np# The following section can be duplicated for building multiple libraries -
#np# make sure to list the correct categories with each library.
#np# replace new_package with the name of the package being autotool'ed here
# ------------------------------------------------------------------------
# Shards library specifications
# ------------------------------------------------------------------------
#np# replace new_package with the name of the package being autotool'ed here
SHARDS_LIB = libshards.a

#np# replace new_package with the name of the package being autotool'ed here
libshards_a_SOURCES = $(SHARDS_CPP)

#np# replace new_package with the name of the package being autotool'ed here
EXTRA_libshards_a_SOURCES =

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
libshards_a_AR = $(ALTERNATE_AR)
else

#np# replace new_package with the name of the package being autotool'ed here
libshards_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
libshards_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
lib_LIBRARIES = libshards.a

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
include_HEADERS = $(SHARDS_H) $(SHARDS_HPP)

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
nodist_execinclude_HEADERS = Shards_config.hpp

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in
