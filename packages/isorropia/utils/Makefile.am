#@HEADER
# ************************************************************************
#
#               Isorropia: Partitioning and Load Balancing Package
#                 Copyright (2006) Sandia Corporation
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
#
# ************************************************************************
#@HEADER

include $(top_builddir)/Makefile.export.isorropia

CORE = \
	$(srcdir)/ispatest_epetra_utils.cpp \
	$(srcdir)/ispatest_lbeval_utils.cpp \
	$(srcdir)/ispatest_utils.cpp \
	$(srcdir)/ispatest_read_distribution.cpp

#Although these headers do not have to be installed, they do need to be
#included in the distribution tarball.
EXTRA_DIST = \
	$(srcdir)/ispatest_epetra_utils.hpp \
	$(srcdir)/ispatest_lbeval_utils.hpp \
	$(srcdir)/ispatest_utils.hpp \
	$(srcdir)/ispatest_read_distribution.hpp

DEFAULT_INCLUDES = -I$(top_srcdir)/utils $(ISORROPIA_INCLUDES)

# ------------------------------------------------------------------------
# Epetra_test library specifications
# ------------------------------------------------------------------------

ISORROPIA_TEST_UTILS_LIB = libisorropia_test_utils.a


libisorropia_test_utils_a_SOURCES = \
	$(CORE)


#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libisorropia_test_utils_a_AR = $(ALTERNATE_AR)
else

libisorropia_test_utils_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# On alpha-dec machines, the extra .o's are needed for Epetra only.
# ------------------------------------------------------------------------

libisorropia_test_utils_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

lib_LIBRARIES = $(ISORROPIA_TEST_UTILS_LIB)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

