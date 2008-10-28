# @HEADER
# ************************************************************************
# 
#                          ThreadPool Package
#                 Copyright (2008) Sandia Corporation
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
# Questions? Contact Carter Edwards (hcedwar@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package ThreadPool subdirectory src built successfully."
	@echo ""

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#np# Make sure to list all source files in one of the following categories.

CORE = \
	$(srcdir)/TPI_pthreads.c \
	$(srcdir)/TPI_Concurrency.c

CORE_H = \
	$(srcdir)/../include/TPI.h \
	$(srcdir)/../include/TPI.hpp

UTIL = \
	$(srcdir)/TPI_Walltime.c


# ------------------------------------------------------------------------
# ThreadPool library specifications
# ------------------------------------------------------------------------
#np# replace new_package with the name of the package being autotool'ed here
THREADPOOL_LIB = libtpi.a

#np# replace new_package with the name of the package being autotool'ed here
THREADPOOL_H = \
	$(CORE_H)

#np# replace new_package with the name of the package being autotool'ed here
libtpi_a_SOURCES = \
	$(CORE) \
	$(UTIL)

#np# replace new_package with the name of the package being autotool'ed here
#EXTRA_libtpi_a_SOURCES =

include $(top_builddir)/Makefile.export.threadpool

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(THREADPOOL_INCLUDES))
else
EXPORT_INCLUDES = $(THREADPOOL_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

if USE_ALTERNATE_AR

libtpi_a_AR = $(ALTERNATE_AR)
else

libtpi_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libtpi_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(THREADPOOL_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
include_HEADERS = $(THREADPOOL_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: New_Package_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(includedir)
#np# replace new_package with the name of the package being autotool'ed here
nodist_execinclude_HEADERS = ../include/ThreadPool_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






