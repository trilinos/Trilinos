# @HEADER
# ************************************************************************
# 
#                Pliris: Parallel Dense Solver Package 
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

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package pliris subdirectory src built successfully."
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
EXTRA_CORE_H = \
        $(srcdir)/block.h \
        $(srcdir)/distribute.h \
        $(srcdir)/exchange_pivots.h \
        $(srcdir)/factor.h \
        $(srcdir)/forwardall.h \
        $(srcdir)/forward.h \
        $(srcdir)/pcomm.h \
        $(srcdir)/perm1.h \
        $(srcdir)/permute.h \
        $(srcdir)/permute_mat.h \
        $(srcdir)/permute_rhs.h \
        $(srcdir)/solve.h \
        $(srcdir)/vars.h \
        $(srcdir)/x_factor.h \
        $(srcdir)/x_solve.h \
        $(srcdir)/xlu_solve.h
#        $(srcdir)/init.h

EXTRA_CORE = \
	$(srcdir)/distribute.c \
        $(srcdir)/exchange_pivots.c \
        $(srcdir)/factor_new.c \
        $(srcdir)/forwardall.c \
        $(srcdir)/forward.c \
        $(srcdir)/pcomm.c \
        $(srcdir)/permute.c \
        $(srcdir)/permute_mat.c \
        $(srcdir)/permute_rhs.c \
        $(srcdir)/solve.c \
        $(srcdir)/x_factor.c \
        $(srcdir)/xlu_solve_new.c \
        $(srcdir)/x_solve.c
#        $(srcdir)/init.c

EXTRA_UTIL_H = \
        $(srcdir)/BLAS_prototypes.h \
        $(srcdir)/cblassp.h \
        $(srcdir)/dblassp.h \
        $(srcdir)/defines.h \
        $(srcdir)/extern_vars.h \
        $(srcdir)/macros.h \
        $(srcdir)/mytime.h \
        $(srcdir)/sblassp.h \
        $(srcdir)/zblassp.h
#        $(srcdir)/clean_code.h

EXTRA_UTIL =  \
	$(srcdir)/clean_code.c \
        $(srcdir)/time.c
#        $(srcdir)/malloc.c \
#        $(srcdir)/my_srand48.c

CORE_H = \
	$(srcdir)/Pliris_ConfigDefs.h

UTIL = \
	$(srcdir)/PlirisVersion.cpp

UTIL_H = \
	$(srcdir)/PlirisVersion.h

EXTRA_EPETRA_COMPAT = \
	$(srcdir)/Pliris.cpp

EXTRA_EPETRA_COMPAT_H = \
	$(srcdir)/Pliris.h

if HAVE_EPETRA

EPETRA_COMPAT = \
	$(EXTRA_EPETRA_COMPAT)

EPETRA_COMPAT_H = \
	$(EXTRA_EPETRA_COMPAT_H)

endif

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
	$(EXTRA_UTIL_H) \
	$(EXTRA_CORE_H)
MPI = \
	$(EXTRA_UTIL) \
        $(EXTRA_CORE)

MPI_EPETRA_H = \
        $(EPETRA_COMPAT_H)
MPI_EPETRA = \
        $(EPETRA_COMPAT)

else

MPI_H =
MPI =

MPI_EPETRA_H =
MPI_EPETRA =

endif

LANG_WRAPS_H = 

LANG_WRAPS = 

DIST_H = 

DIST = 

COMP_H =

COMP = 


#np# List locations of include files within Trilinos but outside of this 
#np# package
#  
include $(top_builddir)/Makefile.export.pliris

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(PLIRIS_INCLUDES))
else
EXPORT_INCLUDES = $(PLIRIS_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# ------------------------------------------------------------------------
# Pliris library specifications
# ------------------------------------------------------------------------

#It would make more sense to say PLIRIS_LIBS here, but that name is
#reserved for the Makefile.export system.
PLIRIS_LIB = libdpliris.a libspliris.a libzpliris.a libcpliris.a

PLIRIS_H = \
	$(CORE_H) \
	$(UTIL_H) \
	$(SERIAL_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H) \
        $(MPI_EPETRA_H)

libdpliris_a_SOURCES = \
	$(CORE) \
	$(UTIL) \
	$(SERIAL) \
	$(DIST) \
	$(COMP) \
	$(MPI) \
        $(MPI_EPETRA)

libspliris_a_SOURCES = \
        $(CORE) \
        $(UTIL) \
        $(SERIAL) \
        $(DIST) \
        $(COMP) \
        $(MPI)

libzpliris_a_SOURCES = \
        $(CORE) \
        $(UTIL) \
        $(SERIAL) \
        $(DIST) \
        $(COMP) \
        $(MPI)

libcpliris_a_SOURCES = \
        $(CORE) \
        $(UTIL) \
        $(SERIAL) \
        $(DIST) \
        $(COMP) \
        $(MPI)

#We only need to list files under the "EXTRA" heading once to
#make sure the files are included in the distribution tarball.
EXTRA_libdpliris_a_SOURCES = \
	$(EXTRA_CORE) \
	$(EXTRA_CORE_H) \
	$(EXTRA_UTIL) \
	$(EXTRA_UTIL_H) \
	$(EXTRA_EPETRA_COMPAT) \
	$(EXTRA_EPETRA_COMPAT_H)

#libdpliris_a_AM_CFLAGS =-DDREAL
libdpliris_a_CFLAGS =-DDREAL
libdpliris_a_CXXFLAGS =-DDREAL

#libspliris_a_AM_CFLAGS =-DSREAL
libspliris_a_CFLAGS =-DSREAL
libspliris_a_CXXFLAGS =-DSREAL

#libcpliris_a_AM_CFLAGS =-DSCPLX
libcpliris_a_CFLAGS =-DSCPLX
libcpliris_a_CXXFLAGS =-DSCPLX

#libzpliris_a_AM_CFLAGS =-DZCPLX
libzpliris_a_CFLAGS =-DZCPLX
libzpliris_a_CXXFLAGS =-DZCPLX

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
libdpliris_a_AR = $(ALTERNATE_AR)
libspliris_a_AR = $(ALTERNATE_AR)
libcpliris_a_AR = $(ALTERNATE_AR)
libzpliris_a_AR = $(ALTERNATE_AR)

else

#np# replace new_package with the name of the package being autotool'ed here
libdpliris_a_AR = $(AR) cru
libspliris_a_AR = $(AR) cru
libcpliris_a_AR = $(AR) cru
libzpliris_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
libdpliris_a_LIBADD = $(XTRALDADD)
libspliris_a_LIBADD = $(XTRALDADD)
libcpliris_a_LIBADD = $(XTRALDADD)
libzpliris_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
lib_LIBRARIES = $(PLIRIS_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
include_HEADERS = $(PLIRIS_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

#np# replace new_package with the name of the package being autotool'ed here
# SPECIAL NOTE: New_Package_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(includedir)
#np# replace new_package with the name of the package being autotool'ed here
nodist_execinclude_HEADERS = Pliris_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

