# @HEADER
# ************************************************************************
# 
#              Finite Element Interface to Linear Solvers
#                 Copyright (2001) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# Questions? Contact Alan Williams (william@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package fei subdirectory support-PETSc built successfully."
	@echo ""

PETSC_INCLUDE=@PETSC_INCLUDE@
PETSC_LIBRARY=@PETSC_LIBRARY@

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside FEI. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

DEFS =
DEFAULT_INCLUDES = -I$(top_srcdir)/interface -I$(top_srcdir) -I$(top_builddir) $(PETSC_INCLUDE)

CORE_H = \
	$(top_srcdir)/support-PETSc/cfei_petsc.h \
	$(top_srcdir)/support-PETSc/fei_PETSc_LinSysCore.hpp

CORE = \
	$(srcdir)/cfei_petsc.cpp \
	$(srcdir)/fei_PETSc_LinSysCore.cpp

# ------------------------------------------------------------------------
# FEI library specifications
# ------------------------------------------------------------------------

FEI_LIB = libfei_petsc.a


FEI_H = \
	$(CORE_H)

libfei_petsc_a_SOURCES = \
	$(CORE) 

EXTRA_libfei_petsc_a_SOURCES = 

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR if USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libfei_petsc_a_AR = $(ALTERNATE_AR)
else

libfei_petsc_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# On alpha-dec machines, the extra .o's are needed for Epetra only.
# ------------------------------------------------------------------------

libfei_petsc_a_LIBADD = $(FEI_TEMPLATE_OBJS) $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(FEI_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(FEI_H) 

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

