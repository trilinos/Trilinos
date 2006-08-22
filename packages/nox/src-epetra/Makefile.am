## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#            NOX: An Object-Oriented Nonlinear Solver Package
#                 Copyright (2002) Sandia Corporation
# 
#            LOCA: Library of Continuation Algorithms Package
#                 Copyright (2005) Sandia Corporation
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
# Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
# Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package nox subdirectory src-epetra built successfully."
	@echo ""

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside NOX. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

EPETRA_H = \
	$(srcdir)/NOX_Epetra.H \
	$(srcdir)/NOX_Epetra_Vector.H \
	$(srcdir)/NOX_Epetra_VectorSpace.H \
	$(srcdir)/NOX_Epetra_VectorSpace_L2.H \
	$(srcdir)/NOX_Epetra_VectorSpace_ScaledL2.H \
	$(srcdir)/NOX_Epetra_MultiVector.H \
	$(srcdir)/NOX_Epetra_Group.H \
	$(srcdir)/NOX_Epetra_Interface_Required.H \
	$(srcdir)/NOX_Epetra_Interface_Jacobian.H \
	$(srcdir)/NOX_Epetra_Interface_Preconditioner.H \
	$(srcdir)/NOX_Epetra_LinearSystem.H \
	$(srcdir)/NOX_Epetra_LinearSystem_AztecOO.H \
        $(srcdir)/NOX_Epetra_MatrixFree.H \
        $(srcdir)/NOX_Epetra_FiniteDifference.H \
        $(srcdir)/NOX_Epetra_FiniteDifferenceColoring.H \
	$(srcdir)/NOX_Epetra_Scaling.H \
	$(srcdir)/NOX_Epetra_BroydenOperator.H

EPETRA = \
	$(srcdir)/NOX_Epetra_Vector.C \
	$(srcdir)/NOX_Epetra_VectorSpace_L2.C \
	$(srcdir)/NOX_Epetra_VectorSpace_ScaledL2.C \
	$(srcdir)/NOX_Epetra_MultiVector.C \
	$(srcdir)/NOX_Epetra_Group.C \
	$(srcdir)/NOX_Epetra_LinearSystem_AztecOO.C \
        $(srcdir)/NOX_Epetra_MatrixFree.C \
        $(srcdir)/NOX_Epetra_FiniteDifference.C \
        $(srcdir)/NOX_Epetra_FiniteDifferenceColoring.C \
	$(srcdir)/NOX_Epetra_Scaling.C \
	$(srcdir)/NOX_Epetra_BroydenOperator.C

if BUILD_NOX_EPETRAEXT_SUPPORT
MODELEVAL = \
        $(srcdir)/NOX_Epetra_ModelEvaluatorInterface.C

MODELEVAL_H = \
        $(srcdir)/NOX_Epetra_ModelEvaluatorInterface.H
endif

# ------------------------------------------------------------------------
# Noxepetra library specifications
# ------------------------------------------------------------------------

libnoxepetra_a_SOURCES = $(EPETRA) $(EPETRA_H) $(MODELEVAL) $(MODELEVAL_H)

# ------------------------------------------------------------------------
# Includes for required and optional library dependencies
# ------------------------------------------------------------------------
include $(top_builddir)/Makefile.export.nox

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(NOX_INCLUDES))
else
EXPORT_INCLUDES = $(NOX_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# ------------------------------------------------------------------------
# When to turn it on
# ------------------------------------------------------------------------

if BUILD_LIBNOXEPETRA

# List of all libraries to install in $(libexecdir)

lib_LIBRARIES = libnoxepetra.a

# List of all headers to install in $(includedir)

include_HEADERS = $(EPETRA_H) $(MODELEVAL_H)

endif

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libnoxepetra_a_AR = $(ALTERNATE_AR)

else

libnoxepetra_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libnoxepetra_a_LIBADD = $(CXXREPOS) $(XTRALDADD)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Files to be deleted by 'make clean'
# ------------------------------------------------------------------------

MOSTLYCLEANFILES = $(CXXREPOS)






