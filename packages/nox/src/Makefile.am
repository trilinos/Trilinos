## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#            NOX: An Object-Oriented Nonlinear Solver Package
#                 Copyright (2002) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#   
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#   
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# 
# Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
# (rppawlo@sandia.gov).
# 
# ************************************************************************
# @HEADER


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside NOX. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------
ABSTRACT_H = \
	$(srcdir)/NOX_Abstract_Vector.H  \
	$(srcdir)/NOX_Abstract_Group.H

LINESEARCH_H = \
	$(srcdir)/NOX_Linesearch_Generic.H \
	$(srcdir)/NOX_Linesearch_FullStep.H \
	$(srcdir)/NOX_Linesearch_Backtrack.H \
	$(srcdir)/NOX_Linesearch_Manager.H \
	$(srcdir)/NOX_Linesearch_MoreThuente.H \
	$(srcdir)/NOX_Linesearch_Polynomial.H 

LINESEARCH = \
	$(srcdir)/NOX_Linesearch_FullStep.C \
	$(srcdir)/NOX_Linesearch_Backtrack.C \
	$(srcdir)/NOX_Linesearch_Manager.C \
	$(srcdir)/NOX_Linesearch_MoreThuente.C \
	$(srcdir)/NOX_Linesearch_Polynomial.C 

DIRECTION_H = \
	$(srcdir)/NOX_Direction_Generic.H \
	$(srcdir)/NOX_Direction_Manager.H \
	$(srcdir)/NOX_Direction_Newton.H \
	$(srcdir)/NOX_Direction_SteepestDescent.H

DIRECTION = \
	$(srcdir)/NOX_Direction_Manager.C \
	$(srcdir)/NOX_Direction_Newton.C \
	$(srcdir)/NOX_Direction_SteepestDescent.C

PARAMETER_H = \
	$(srcdir)/NOX_Parameter_Entry.H \
	$(srcdir)/NOX_Parameter_List.H 

PARAMETER = \
	$(srcdir)/NOX_Parameter_Entry.C \
	$(srcdir)/NOX_Parameter_List.C

SOLVER_H = \
	$(srcdir)/NOX_Solver_Generic.H \
	$(srcdir)/NOX_Solver_Manager.H \
	$(srcdir)/NOX_Solver_LineSearch.H \
	$(srcdir)/NOX_Solver_TrustRegion.H \
	$(srcdir)/NOX_Solver_NonlinearCG.H 

SOLVER = \
	$(srcdir)/NOX_Solver_Manager.C \
	$(srcdir)/NOX_Solver_LineSearch.C \
	$(srcdir)/NOX_Solver_TrustRegion.C \
	$(srcdir)/NOX_Solver_NonlinearCG.C 

STATUS_H = \
	$(srcdir)/NOX_Status_Test.H \
	$(srcdir)/NOX_Status_AbsResid.H \
	$(srcdir)/NOX_Status_Combo.H \
	$(srcdir)/NOX_Status_MaxResid.H \
	$(srcdir)/NOX_Status_RelResid.H \
	$(srcdir)/NOX_Status_MaxIters.H \
	$(srcdir)/NOX_Status_RelUpdate.H 

STATUS = \
	$(srcdir)/NOX_Status_Test.C \
	$(srcdir)/NOX_Status_AbsResid.C \
	$(srcdir)/NOX_Status_Combo.C \
	$(srcdir)/NOX_Status_MaxResid.C \
	$(srcdir)/NOX_Status_RelResid.C \
	$(srcdir)/NOX_Status_MaxIters.C \
	$(srcdir)/NOX_Status_RelUpdate.C 

UTILS_H = \
	$(srcdir)/NOX_Utils.H 

UTILS = \
	$(srcdir)/NOX_Utils.C


EPETRA_H = \
	$(srcdir)/NOX_Epetra_Vector.H \
	$(srcdir)/NOX_Epetra_Group.H \
	$(srcdir)/NOX_Epetra_SharedJacobian.H \
        $(srcdir)/NOX_Epetra_MatrixFree.H \
        $(srcdir)/NOX_Epetra_FiniteDifference.H \
	$(srcdir)/NOX_Epetra_Interface.H 

EPETRA = \
	$(srcdir)/NOX_Epetra_Vector.C \
	$(srcdir)/NOX_Epetra_Group.C \
	$(srcdir)/NOX_Epetra_SharedJacobian.C \
        $(srcdir)/NOX_Epetra_MatrixFree.C \
        $(srcdir)/NOX_Epetra_FiniteDifference.C

# ------------------------------------------------------------------------
# Nox library specifications
# ------------------------------------------------------------------------

NOX_LIB = libnox.a

NOX_H = \
	$(ABSTRACT_H) \
	$(LINESEARCH_H) \
	$(DIRECTION_H) \
	$(PARAMETER_H) \
	$(SOLVER_H) \
	$(STATUS_H) \
	$(UTILS_H) 


libnox_a_SOURCES = \
	$(LINESEARCH) \
	$(DIRECTION) \
	$(PARAMETER) \
	$(SOLVER) \
	$(STATUS) \
	$(UTILS) 


# ------------------------------------------------------------------------
# Noxepetra library specifications
# ------------------------------------------------------------------------

if BUILD_NOXEPETRA

NOXEPETRA_LIB = libnoxepetra.a

NOXEPETRA_H = $(EPETRA_H)

libnoxepetra_a_SOURCES = $(EPETRA)

libnoxepetra_a_CXXFLAGS = $(TRILINOS_CXXFLAGS)

endif

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libnox_a_AR = $(ALTERNATE_AR)
libnoxepetra_a_AR = $(ALTERNATE_AR)

else

libnox_a_AR = $(AR) cru
libnoxepetra_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libnox_a_LIBADD = $(XTRALDADD)
libnoxepetra_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(NOX_LIB) $(NOXEPETRA_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(pkgincludedir)
# ------------------------------------------------------------------------

pkginclude_HEADERS = $(NOX_H) $(NOXEPETRA_H)

# ------------------------------------------------------------------------
# Special stuff to install in our special $(pkgexecincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: NOX_Config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

pkgexecincludedir = $(exec_prefix)/include
pkgexecinclude_HEADERS = NOX_Config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






