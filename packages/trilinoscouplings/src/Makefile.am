# @HEADER
# ************************************************************************
# 
#                     TrilinosCouplings Package
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
	@echo "Trilinos package trilinoscouplings subdirectory src built successfully."
	@echo ""

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside TrilinosCouplings. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
# There is, at present, no distinction between most of the categories, 
# other than those which start with EXTRA_.  All categories which do not
# start with EXTRA_ are included in all builds.  Categories which start
# with EXTRA_ are included in some builds (typically controlled by 
# --with-xyz or --without-xyz)
#
# Make sure to list all source files in one of the following categories.
CORE_H = \
	$(srcdir)/TrilinosCouplings_ConfigDefs.h

CORE = \
	$(srcdir)/TrilinosCouplings_Version.h 

UTIL_H =

UTIL = 

#
#  Files which are included in select builds are typically have
#  names that start with EXTRA
#

# Because TrilinosCouplings interfaces are only build when all required
# packages are available, interface files must be listed as conditionally
# compiled files (ie with an EXTRA_ prefix).

EXTRA_NOX_ML_H = \
    $(srcdir)/ml/NonlinML/nlnml_ConstrainedMultiLevelOperator.H \
    $(srcdir)/ml/NonlinML/nlnml_finelevelnoxinterface.H         \
    $(srcdir)/ml/NonlinML/nlnml_coarselevelnoxinterface.H       \
    $(srcdir)/ml/NonlinML/nlnml_preconditioner.H                \
    $(srcdir)/ml/NonlinML/nlnml_preconditioner_utils.H          \
    $(srcdir)/ml/NonlinML/nlnml_linearsystem.H                  \
    $(srcdir)/ml/NonlinML/nlnml_nonlinearlevel.H                \
    $(srcdir)/ml/NonlinML/nlnml_prepostoperator.H

EXTRA_NOX_ML = \
    $(srcdir)/ml/NonlinML/nlnml_prepostoperator.cpp                     \
    $(srcdir)/ml/NonlinML/nlnml_ConstrainedMultiLevelOperator.cpp       \
    $(srcdir)/ml/NonlinML/nlnml_linearsystem.cpp                        \
    $(srcdir)/ml/NonlinML/nlnml_nonlinearlevel.cpp                      \
    $(srcdir)/ml/NonlinML/nlnml_coarselevelnoxinterface.cpp             \
    $(srcdir)/ml/NonlinML/nlnml_preconditioner1.cpp                     \
    $(srcdir)/ml/NonlinML/nlnml_preconditioner_utils.cpp

if BUILD_NOX_ML

NOX_ML_H = \
	$(EXTRA_NOX_ML_H)
NOX_ML = \
	$(EXTRA_NOX_ML)

else

NOX_ML_H =
NOX_ML = 

endif

EXTRA_IFPACK_ML_H = \
	$(srcdir)/ifpack/Ifpack_ML.h 

EXTRA_IFPACK_ML = \
	$(srcdir)/ifpack/Ifpack_ML.cpp

if BUILD_IFPACK_ML 

IFPACK_ML_H = \
	$(EXTRA_IFPACK_ML_H)
IFPACK_ML = \
	$(EXTRA_IFPACK_ML)

else

IFPACK_ML_H =
IFPACK_ML =

endif

EXTRA_EPETRAEXT_ISORROPIA_H = \
	$(srcdir)/epetraext/EpetraExt_Isorropia_CrsGraph.h

EXTRA_EPETRAEXT_ISORROPIA = \
	$(srcdir)/epetraext/EpetraExt_Isorropia_CrsGraph.cpp

if BUILD_EPETRAEXT_ISORROPIA 

EPETRAEXT_ISORROPIA_H = \
	$(EXTRA_EPETRAEXT_ISORROPIA_H)
EPETRAEXT_ISORROPIA = \
	$(EXTRA_EPETRAEXT_ISORROPIA)

else      

EPETRAEXT_ISORROPIA_H =
EPETRAEXT_ISORROPIA = 

endif

EXTRA_EPETRAEXT_AMESOS_H = \
	$(srcdir)/epetraext/EpetraExt_AmesosBTF_CrsGraph.h \
	$(srcdir)/epetraext/EpetraExt_AmesosBTF_CrsMatrix.h 

EXTRA_EPETRAEXT_AMESOS = \
	$(srcdir)/epetraext/EpetraExt_AmesosBTF_CrsGraph.cpp \
	$(srcdir)/epetraext/EpetraExt_AmesosBTF_CrsMatrix.cpp

if BUILD_EPETRAEXT_AMESOS

EPETRAEXT_AMESOS_H = \
	$(EXTRA_EPETRAEXT_AMESOS_H)

EPETRAEXT_AMESOS = \
	$(EXTRA_EPETRAEXT_AMESOS)

else

EPETRAEXT_AMESOS_H =
EPETRAEXT_AMESOS =

endif

EXTRA_EPETRAEXT_AMESOS_ISORROPIA_H = \
	$(srcdir)/epetraext/EpetraExt_AmesosBTFGlobal_LinearProblem.h

EXTRA_EPETRAEXT_AMESOS_ISORROPIA = \
	$(srcdir)/epetraext/EpetraExt_AmesosBTFGlobal_LinearProblem.cpp

if BUILD_EPETRAEXT_AMESOS_ISORROPIA

EPETRAEXT_AMESOS_ISORROPIA_H = \
	$(EXTRA_EPETRAEXT_AMESOS_ISORROPIA_H)

EPETRAEXT_AMESOS_ISORROPIA = \
	$(EXTRA_EPETRAEXT_AMESOS_ISORROPIA)

else

EPETRAEXT_AMESOS_ISORROPIA_H =
EPETRAEXT_AMESOS_ISORROPIA =

endif

##conditional NOX dependencies
#EXTRA_NOX_ML_H = \
#   -I$(top_srcdir)/../nox/src -I$(top_srcdir)/../nox/src-epetra -I$(top_builddir)/../nox/src
#if BUILD_NOX_ML
#NOX_ML_H = \
#  $(EXTRA_NOX_ML_H)
#else
#NOX_ML_H =
#endif

#AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src

# ------------------------------------------------------------------------
# TrilinosCouplings library specifications
# ------------------------------------------------------------------------
#np# replace trilinoscouplings with the name of the package being autotool'ed here
TRILINOSCOUPLINGS_LIB = libtrilinoscouplings.a

#np# replace trilinoscouplings with the name of the package being autotool'ed here
TRILINOSCOUPLINGS_H = \
	$(CORE_H) \
	$(UTIL_H) \
	$(IFPACK_ML_H) \
	$(NOX_ML_H) \
	$(EPETRAEXT_ISORROPIA_H) \
	$(EPETRAEXT_AMESOS_H) \
	$(EPETRAEXT_AMESOS_ISORROPIA_H)      

#np# replace trilinoscouplings with the name of the package being autotool'ed here
libtrilinoscouplings_a_SOURCES = \
	$(CORE) \
	$(UTIL) \
	$(IFPACK_ML) \
	$(NOX_ML) \
	$(EPETRAEXT_ISORROPIA) \
	$(EPETRAEXT_AMESOS) \
	$(EPETRAEXT_AMESOS_ISORROPIA)


EXTRA_libtrilinoscouplings_a_SOURCES = \
	$(EXTRA_IFPACK_ML) \
	$(EXTRA_IFPACK_ML_H) \
	$(EXTRA_NOX_ML) \
	$(EXTRA_NOX_ML_H) \
	$(EXTRA_EPETRAEXT_ISORROPIA) \
	$(EXTRA_EPETRAEXT_ISORROPIA_H) \
	$(EXTRA_EPETRAEXT_AMESOS) \
	$(EXTRA_EPETRAEXT_AMESOS_H) \
	$(EXTRA_EPETRAEXT_AMESOS_ISORROPIA) \
	$(EXTRA_EPETRAEXT_AMESOS_ISORROPIA_H)

include $(top_builddir)/Makefile.export.trilinoscouplings

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(TRILINOSCOUPLINGS_INCLUDES))
else
EXPORT_INCLUDES = $(TRILINOSCOUPLINGS_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)


# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libtrilinoscouplings_a_AR = $(ALTERNATE_AR)
else

libtrilinoscouplings_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libtrilinoscouplings_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(TRILINOSCOUPLINGS_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(TRILINOSCOUPLINGS_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: TrilinosCouplings_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(includedir)
nodist_execinclude_HEADERS = TrilinosCouplings_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






