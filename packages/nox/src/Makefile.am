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
ABSTRACT = \
	$(srcdir)/NOX_Abstract_Group.C

ABSTRACT_H = \
	$(srcdir)/NOX_Abstract_Vector.H  \
	$(srcdir)/NOX_Abstract_Group.H

LINESEARCH_H = \
	$(srcdir)/NOX_LineSearch_Generic.H \
	$(srcdir)/NOX_LineSearch_FullStep.H \
	$(srcdir)/NOX_LineSearch_Backtrack.H \
	$(srcdir)/NOX_LineSearch_Manager.H \
	$(srcdir)/NOX_LineSearch_MoreThuente.H \
	$(srcdir)/NOX_LineSearch_Polynomial.H \
	$(srcdir)/NOX_LineSearch_Secant.H 

LINESEARCH = \
	$(srcdir)/NOX_LineSearch_FullStep.C \
	$(srcdir)/NOX_LineSearch_Backtrack.C \
	$(srcdir)/NOX_LineSearch_Manager.C \
	$(srcdir)/NOX_LineSearch_MoreThuente.C \
	$(srcdir)/NOX_LineSearch_Polynomial.C \
	$(srcdir)/NOX_LineSearch_Secant.C 

DIRECTION_H = \
	$(srcdir)/NOX_Direction_Generic.H \
	$(srcdir)/NOX_Direction_Manager.H \
	$(srcdir)/NOX_Direction_Newton.H \
	$(srcdir)/NOX_Direction_NonlinearCG.H \
	$(srcdir)/NOX_Direction_SteepestDescent.H

DIRECTION = \
	$(srcdir)/NOX_Direction_Manager.C \
	$(srcdir)/NOX_Direction_Newton.C \
	$(srcdir)/NOX_Direction_NonlinearCG.C \
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
	$(srcdir)/NOX_Solver_LineSearchBased.H \
	$(srcdir)/NOX_Solver_TrustRegionBased.H

SOLVER = \
	$(srcdir)/NOX_Solver_Manager.C \
	$(srcdir)/NOX_Solver_LineSearchBased.C \
	$(srcdir)/NOX_Solver_TrustRegionBased.C

STATUSTEST_H = \
	$(srcdir)/NOX_StatusTest_Generic.H \
	$(srcdir)/NOX_StatusTest_NormF.H \
	$(srcdir)/NOX_StatusTest_NormUpdate.H \
	$(srcdir)/NOX_StatusTest_Combo.H \
	$(srcdir)/NOX_StatusTest_MaxIters.H \
        $(srcdir)/NOX_StatusTest_NormWRMS.H 

STATUSTEST = \
	$(srcdir)/NOX_StatusTest_Generic.C \
	$(srcdir)/NOX_StatusTest_NormF.C \
	$(srcdir)/NOX_StatusTest_NormUpdate.C \
	$(srcdir)/NOX_StatusTest_Combo.C \
	$(srcdir)/NOX_StatusTest_MaxIters.C \
        $(srcdir)/NOX_StatusTest_NormWRMS.C  

UTILS_H = \
	$(srcdir)/NOX_Utils.H 

UTILS = \
	$(srcdir)/NOX_Utils.C

# ------------------------------------------------------------------------
# Nox library specifications
# ------------------------------------------------------------------------

NOX_LIB = libnox.a

NOX_H = \
	$(srcdir)/NOX.H \
	$(srcdir)/NOX_Description.H \
	$(srcdir)/NOX_Common.H \
	$(ABSTRACT_H) \
	$(LINESEARCH_H) \
	$(DIRECTION_H) \
	$(PARAMETER_H) \
	$(SOLVER_H) \
	$(STATUSTEST_H) \
	$(UTILS_H) 


libnox_a_SOURCES = \
	$(ABSTRACT) \
	$(LINESEARCH) \
	$(DIRECTION) \
	$(PARAMETER) \
	$(SOLVER) \
	$(STATUSTEST) \
	$(UTILS) 


# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libnox_a_AR = $(ALTERNATE_AR)

else

libnox_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libnox_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(NOX_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(pkgincludedir)
# ------------------------------------------------------------------------

pkginclude_HEADERS = $(NOX_H)

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






