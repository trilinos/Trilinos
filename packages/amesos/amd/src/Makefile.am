# @HEADER
# ************************************************************************
#
# The copyright to the AMD package belongs to Tim Davis (davis@cse.ufl.edu) 
# This is not the official version of this AMD package (see availability
# below for the official version).  It has been modified as follows for
# use within Amesos, the sparse directs olver package in Trilinos:
#   1)  Only the integer interface is compiled and linked into libamd.a
#   
# 
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
#
# AMD Version 1.0 (Apr. 30, 2003),  Copyright (c) 2003 by Timothy A.
# Davis, Patrick R. Amestoy, and Iain S. Duff.  All Rights Reserved.
# 
# AMD License:
# 
#     Your use or distribution of AMD or any modified version of
#     AMD implies that you agree to this License.
# 
#     THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
#     EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
# 
#     Permission is hereby granted to use or copy this program, provided
#     that the Copyright, this License, and the Availability of the original
#     version is retained on all copies.  User documentation of any code that
#     uses AMD or any modified version of AMD code must cite the
#     Copyright, this License, the Availability note, and "Used by permission."
#     Permission to modify the code and to distribute modified code is granted,
#     provided the Copyright, this License, and the Availability note are
#     retained, and a notice that the code was modified is included.  This
#     software was developed with support from the National Science Foundation,
#     and is provided to you free of charge.
# 
# Availability:
# 
#     http://www.cise.ufl.edu/research/sparse/amd
# 
# ************************************************************************
# @HEADER

# Amesos uses amd, but we don't make amd available outside Amesos_Klu 
# at this time.  
#

# I uncommented the following lines so that the files are included in
# the distribution tarball.  We could probably just list them in an 
# EXTRA_DIST macro, but there isn't yet a precedent for this with headers. JW

 CORE_H = \
	$(srcdir)/amd_internal.h  \
	$(srcdir)/amd.h 

CORE = \
	$(srcdir)/amd_aat.c \
	$(srcdir)/amd_1.c \
	$(srcdir)/amd_2.c \
	$(srcdir)/amd_dump.c  \
	$(srcdir)/amd_postorder.c \
	$(srcdir)/amd_post_tree.c \
	$(srcdir)/amd_defaults.c  \
	$(srcdir)/amd_order.c  \
	$(srcdir)/amd_control.c \
	$(srcdir)/amd_info.c \
	$(srcdir)/amd_valid.c  

AM_CPPFLAGS = -DDINT 

AMD_LIB = libamd.a

AMD_H = \
	$(CORE_H)

#np# replace amd with the name of the package being autotool'ed here
libamd_a_SOURCES = \
	$(CORE) 

#np# replace amd with the name of the package being autotool'ed here
EXTRA_libamd_a_SOURCES =


#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

#np# replace amd with the name of the package being autotool'ed here
libamd_a_AR = $(ALTERNATE_AR)
else

#np# replace amd with the name of the package being autotool'ed here
libamd_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

#np# replace amd with the name of the package being autotool'ed here
libamd_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

#np# replace amd with the name of the package being autotool'ed here
lib_LIBRARIES = $(AMD_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

#np# replace amd with the name of the package being autotool'ed here
include_HEADERS = $(AMD_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------


# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = 

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






