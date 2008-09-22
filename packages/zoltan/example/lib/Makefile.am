# @HEADER
##############################################################################
# Zoltan Library for Parallel Applications                                   #
# Copyright (c) 2000,2001,2002, Sandia National Laboratories.                #
# This software is distributed under the GNU Lesser General Public License.  #
# For more info, see the README file in the top-level Zoltan directory.      #
##############################################################################
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:

UTIL_H = \
	$(srcdir)/exzoltan.h

UTIL = \
	$(srcdir)/exmesh.c \
	$(srcdir)/exmisc.c \
	$(srcdir)/exphg.c

# ------------------------------------------------------------------------
# ExZoltan library specifications
# ------------------------------------------------------------------------
EXZOLTAN_LIB = libexzoltan.a

EXZOLTAN_H = \
	$(UTIL_H)

libexzoltan_a_SOURCES = \
	$(UTIL)

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR
libexzoltan_a_AR = $(ALTERNATE_AR)
else
libexzoltan_a_AR = $(AR) cru
endif

include $(top_builddir)/Makefile.export.zoltan

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(ZOLTAN_INCLUDES))
else
EXPORT_INCLUDES = $(ZOLTAN_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libexzoltan_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(EXZOLTAN_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(EXZOLTAN_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include
#np# replace new_package with the name of the package being autotool'ed here
#nodist_execinclude_HEADERS = Zoltan_config.h

#SUBDIRS = 

EXTRA_DIST = makefile
# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

