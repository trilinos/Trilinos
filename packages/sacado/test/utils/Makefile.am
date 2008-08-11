## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#                           Sacado Package
#                 Copyright (2006) Sandia Corporation
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
# Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
# (etphipp@sandia.gov).
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside Sacado. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

UTILS_H = \
	$(srcdir)/Sacado_Random.hpp 

UTILS = \
	$(srcdir)/Sacado_Random.cpp

# ------------------------------------------------------------------------
# Utils library specifications
# ------------------------------------------------------------------------
UTILS_LIB = libsacadoutils.a

libsacadoutils_a_SOURCES = $(UTILS) $(UTILS_H)

libsacadoutils_a_CXXFLAGS = \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libsacadoutils_a_AR = $(ALTERNATE_AR)

else

libsacadoutils_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libsacadoutils_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries not to be installed
# ------------------------------------------------------------------------

noinst_LIBRARIES = $(UTILS_LIB)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 
