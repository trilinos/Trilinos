# @HEADER
# ************************************************************************
# 
#     EpetraExt: Epetra Extended - Linear Algebra Services Package 
#                 Copyright (2001) Sandia Corporation
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

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside EpetraExt. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

TRANSFORM_H = \
	$(srcdir)/Epetra_Transform.h \
        $(srcdir)/Epetra_Transform_Composite.h \
	$(srcdir)/EDT_CrsGraph_View.h \
	$(srcdir)/EDT_CrsMatrix_View.h \
	$(srcdir)/EDT_MultiVector_View.h \
	$(srcdir)/EDT_CrsGraph_Transpose.h \
	$(srcdir)/EDT_RowMatrix_Transpose.h \
        $(srcdir)/EDT_CrsGraph_MapColoring.h \
        $(srcdir)/EDT_CrsGraph_MapColoringIndex.h \
	$(srcdir)/EDT_LinearProblem_StaticCondensation.h \
	$(srcdir)/EDT_LinearProblem_CrsSingletonFilter.h \
        $(srcdir)/EDT_CrsGraph_SymmRCM.h \
        $(srcdir)/EDT_LinearProblem_GraphTrans.h

TRANSFORM = \
	$(srcdir)/EDT_CrsGraph_View.cpp \
	$(srcdir)/EDT_CrsMatrix_View.cpp \
	$(srcdir)/EDT_MultiVector_View.cpp \
	$(srcdir)/EDT_CrsGraph_Transpose.cpp \
	$(srcdir)/EDT_RowMatrix_Transpose.cpp \
        $(srcdir)/EDT_CrsGraph_MapColoring.cpp \
        $(srcdir)/EDT_CrsGraph_MapColoringIndex.cpp \
	$(srcdir)/EDT_LinearProblem_StaticCondensation.cpp \
	$(srcdir)/EDT_LinearProblem_CrsSingletonFilter.cpp \
        $(srcdir)/EDT_CrsGraph_SymmRCM.cpp \
        $(srcdir)/EDT_LinearProblem_GraphTrans.cpp

BTF_H = \
        $(srcdir)/EDT_CrsGraph_BTF.h

BTF = \
        $(srcdir)/EDT_CrsGraph_BTF.cpp \
        $(srcdir)/btf/concmp.f \
        $(srcdir)/btf/genbtf.f \
        $(srcdir)/btf/corsum.f \
        $(srcdir)/btf/finsum.f \
        $(srcdir)/btf/fnrsum.f \
        $(srcdir)/btf/izero.f \
        $(srcdir)/btf/maxmatch.f \
        $(srcdir)/btf/rectblk.f \
        $(srcdir)/btf/mmc13e.f \
        $(srcdir)/btf/prtivs.f \
        $(srcdir)/btf/strchk.f \
        $(srcdir)/btf/mattrans.f

EXTRA_EPETRAEXT_ZOLTAN_H = \
	$(srcdir)/EDT_CrsGraph_Zoltan.h \
	$(srcdir)/EDT_LinearProblem_Zoltan.h \
	$(srcdir)/EDT_CrsGraph_ZoltanOrder.h \
	$(srcdir)/Epetra_ZoltanQuery.h

EXTRA_EPETRAEXT_ZOLTAN = \
	$(srcdir)/EDT_CrsGraph_Zoltan.cpp \
	$(srcdir)/EDT_LinearProblem_Zoltan.cpp \
	$(srcdir)/EDT_CrsGraph_ZoltanOrder.cpp \
	$(srcdir)/Epetra_ZoltanQuery.cpp

if HAVE_EPETRAEXT_ZOLTAN

EPETRAEXT_ZOLTAN_H = \
        $(EXTRA_EPETRAEXT_ZOLTAN_H)
EPETRAEXT_ZOLTAN = \
        $(EXTRA_EPETRAEXT_ZOLTAN)

else

EPETRAEXT_ZOLTAN_H =
EPETRAEXT_ZOLTAN =

endif

EXTRA_MPI_H = \
	$(srcdir)/EDT_CrsGraph_Overlap.h \
	$(EPETRAEXT_ZOLTAN_H)

EXTRA_MPI = \
	$(srcdir)/EDT_CrsGraph_Overlap.cpp \
	$(EPETRAEXT_ZOLTAN)
if HAVE_MPI

MPI_H = \
        $(EXTRA_MPI_H)
MPI = \
        $(EXTRA_MPI)

else

MPI_H =
MPI =

endif

# ------------------------------------------------------------------------
# EpetraExt Tranform library specifications
# ------------------------------------------------------------------------

libepetraext_trans_a_SOURCES = \
	$(TRANSFORM) \
        $(BTF) \
	$(MPI)

#If these are not listed, the files will not be included in
#the distribution tarball.
EXTRA_libepetraext_trans_a_SOURCES = \
	$(EXTRA_MPI) \
	$(EXTRA_MPI_H) \
	$(EXTRA_EPETRAEXT_ZOLTAN) \
	$(EXTRA_EPETRAEXT_ZOLTAN_H)

if BUILD_TRANSFORM

# libraries to be installed in $(libexecdir)
EPETRAEXT_TRANSFORM_LIB = libepetraext_trans.a
EPETRAEXT_TRANSFORM_H = \
	$(TRANSFORM_H) \
        $(BTF_H) \
	$(MPI_H)

else

EPETRAEXT_TRANSFORM_LIB =
EPETRAEXT_TRANSFORM_H =

endif

#Eventaully the Zoltan code will not be available in Trilinos/external/
#Until then, we will set up this path automatically.  Users who get the
#Zoltan source from somewhere else will need to specify the path.
if HAVE_EPETRAEXT_ZOLTAN

EPETRAEXT_ZOLTAN_INCDIR=-I$(top_srcdir)/../../external/Zoltan_C++

else

EPETRAEXT_ZOLTAN_INCDIR=

endif

#We need the following to support separate build directories
AM_CPPFLAGS=-I$(top_builddir)/../epetra/src -I$(top_srcdir)/../epetra/src \
            $(EPETRAEXT_ZOLTAN_INCDIR)

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libepetraext_trans_a_AR = $(ALTERNATE_AR)

else

libepetraext_trans_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# On alpha-dec machines, the extra .o's are needed for Epetra only.
# ------------------------------------------------------------------------

libepetraext_trans_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(EPETRAEXT_TRANSFORM_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(EPETRAEXT_TRANSFORM_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Epetra_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS =../EpetraExt_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Files to be deleted by 'make clean'
# ------------------------------------------------------------------------

MOSTLYCLEANFILES = $(TRILINOS_TEMPLATE_OBJS)

