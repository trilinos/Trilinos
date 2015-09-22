#@HEADER
# ************************************************************************
# 
#                     Pamgen Package
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
# Questions? Contact David Hensinger(dmhensi@sandia.gov) 
# 
# ************************************************************************
#@HEADER

include $(top_builddir)/Makefile.export.pamgen

# The following line helps the test harness recover from build errors.

all-local: 

EXEEXT = .exe

noinst_PROGRAMS = pamgen_test_driver 

pamgen_test_driver_SOURCES = $(srcdir)/pamgen_test_driver.c $(srcdir)/getopts.C

# EXTRA_DIST will force these files to be included in the distribution tarball
EXTRA_DIST = $(srcdir)/getopts.h $(srcdir)/copyright.h

#np# Specify locations of include files within Trilinos, but outside np
#np# Add packages that np is dependent on; remove epetra if np does not need it

AM_CPPFLAGS = $(PAMGEN_INCLUDES)

#np# - List libraries that np is dependent on.  The DEPENDENCIES primary is
#np# optional.
pamgen_test_driver_DEPENDENCIES = \
	$(top_builddir)/src/libpamgen.a 

#np# List libraries that np is dependent on.
pamgen_test_driver_LDADD = $(PAMGEN_LIBS)

#np# Set "SUBDIRS= (list all subdirectories that are autotool'ed here)"
#np# List all subdirectories even if one or more are configured/built only
#np# conditionally.
#np# For packages whose examples are located in this directory, rather than
#np# in a subdirectory of this directory, use the Makefile.am located in
#np# example1/ rather than this file.
# SUBDIRS = example1
SUBDIRS = 
