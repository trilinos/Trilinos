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
# Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
# (rppawlo@sandia.gov), Sandia National Laboratories.
# 
# ************************************************************************
# @HEADER


# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = libnox.a

# ------------------------------------------------------------------------
# Nox library specifications
# ------------------------------------------------------------------------


libnox_a_SOURCES = \
	$(ABSTRACT) \
	$(LINESEARCH) \
	$(DIRECTION) \
	$(PARAMETER) \
	$(SOLVER) \
	$(STATUSTEST) \
	$(UTILS) \
	$(MULTIVECS)

libnox_a_CXXFLAGS = -I$(top_srcdir)/src $(TRILINOS_CXXFLAGS)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = \
	$(srcdir)/NOX.H \
	$(srcdir)/NOX_Common.H \
	$(ABSTRACT_H) \
	$(LINESEARCH_H) \
	$(DIRECTION_H) \
	$(PARAMETER_H) \
	$(SOLVER_H) \
	$(STATUSTEST_H) \
	$(UTILS_H) \
	$(MULTIVECS_H)

# ------------------------------------------------------------------------
# Parts of the code that are not ready to be included in nightly release
# ------------------------------------------------------------------------


if BUILD_PRERELEASE

nodist_include_HEADERS = \
	$(srcdir)/NOX_Direction_Tensor.H \
	$(srcdir)/NOX_Solver_InexactTrustRegionBased.H \
	$(srcdir)/NOX_Solver_TensorBased.H \
	$(srcdir)/NOX_Solver_TensorBasedTest.H \
	$(srcdir)/NOX_LineSearch_NonlinearCG.H \
	$(srcdir)/NOX_Direction_NonlinearCG.H \
	$(srcdir)/NOX_Direction_QuasiNewton.H \
	$(srcdir)/NOX_Direction_Broyden.H \
	$(srcdir)/NOX_Direction_ModifiedNewton.H \
        $(srcdir)/NOX_Direction_Utils_InexactNewton.H


nodist_libnox_a_SOURCES = \
	$(srcdir)/NOX_Direction_Tensor.C \
	$(srcdir)/NOX_Solver_InexactTrustRegionBased.C \
	$(srcdir)/NOX_Solver_TensorBased.C \
	$(srcdir)/NOX_Solver_TensorBasedTest.C \
	$(srcdir)/NOX_LineSearch_NonlinearCG.C \
	$(srcdir)/NOX_Direction_NonlinearCG.C \
	$(srcdir)/NOX_Direction_QuasiNewton.C \
	$(srcdir)/NOX_Direction_Broyden.C \
	$(srcdir)/NOX_Direction_ModifiedNewton.C \
        $(srcdir)/NOX_Direction_Utils_InexactNewton.C

endif


# ------------------------------------------------------------------------
# Extra stuff to be included in the distribution
# ------------------------------------------------------------------------

EXTRA_DIST = \
	$(srcdir)/NOX_Description.H 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(pkgexecincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: NOX_Config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include/
nodist_execinclude_HEADERS = NOX_Config.h

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside NOX. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

ABSTRACT = \
	$(srcdir)/NOX_Abstract_Vector.C \
	$(srcdir)/NOX_Abstract_Group.C

ABSTRACT_H = \
	$(srcdir)/NOX_Abstract_Vector.H \
	$(srcdir)/NOX_Abstract_Group.H

LINESEARCH_H = \
	$(srcdir)/NOX_LineSearch_Generic.H \
	$(srcdir)/NOX_LineSearch_Utils_Slope.H \
	$(srcdir)/NOX_LineSearch_Utils_Printing.H \
	$(srcdir)/NOX_LineSearch_Utils_Counters.H \
	$(srcdir)/NOX_LineSearch_FullStep.H \
	$(srcdir)/NOX_LineSearch_Backtrack.H \
	$(srcdir)/NOX_LineSearch_Manager.H \
	$(srcdir)/NOX_LineSearch_MoreThuente.H \
	$(srcdir)/NOX_LineSearch_Polynomial.H 

LINESEARCH = \
	$(srcdir)/NOX_LineSearch_Utils_Slope.C \
	$(srcdir)/NOX_LineSearch_Utils_Printing.C \
	$(srcdir)/NOX_LineSearch_Utils_Counters.C \
	$(srcdir)/NOX_LineSearch_FullStep.C \
	$(srcdir)/NOX_LineSearch_Backtrack.C \
	$(srcdir)/NOX_LineSearch_Manager.C \
	$(srcdir)/NOX_LineSearch_MoreThuente.C \
	$(srcdir)/NOX_LineSearch_Polynomial.C 

DIRECTION_H = \
	$(srcdir)/NOX_Direction_Generic.H \
	$(srcdir)/NOX_Direction_Manager.H \
	$(srcdir)/NOX_Direction_Newton.H \
	$(srcdir)/NOX_Direction_SteepestDescent.H 

DIRECTION = \
	$(srcdir)/NOX_Direction_Manager.C \
	$(srcdir)/NOX_Direction_Generic.C \
	$(srcdir)/NOX_Direction_Newton.C \
	$(srcdir)/NOX_Direction_SteepestDescent.C 

PARAMETER_H = \
	$(srcdir)/NOX_Parameter_Entry.H \
        $(srcdir)/NOX_Parameter_Arbitrary.H \
        $(srcdir)/NOX_Parameter_UserNorm.H \
        $(srcdir)/NOX_Parameter_MeritFunction.H \
        $(srcdir)/NOX_Parameter_PrePostOperator.H \
	$(srcdir)/NOX_Parameter_DirectionConstructor.H \
	$(srcdir)/NOX_Parameter_DirectionConstructorT.H \
	$(srcdir)/NOX_Parameter_LineSearchConstructor.H \
	$(srcdir)/NOX_Parameter_LineSearchConstructorT.H \
	$(srcdir)/NOX_Parameter_AnyPtr.H \
	$(srcdir)/NOX_Parameter_List.H 

PARAMETER = \
	$(srcdir)/NOX_Parameter_Entry.C \
        $(srcdir)/NOX_Parameter_PrePostOperator.C \
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
	$(srcdir)/NOX_StatusTest_FiniteValue.H \
	$(srcdir)/NOX_StatusTest_NormUpdate.H \
	$(srcdir)/NOX_StatusTest_Combo.H \
	$(srcdir)/NOX_StatusTest_MaxIters.H \
        $(srcdir)/NOX_StatusTest_NormWRMS.H \
        $(srcdir)/NOX_StatusTest_Stagnation.H \
        $(srcdir)/NOX_StatusTest_Divergence.H 

STATUSTEST = \
	$(srcdir)/NOX_StatusTest_Generic.C \
	$(srcdir)/NOX_StatusTest_NormF.C \
	$(srcdir)/NOX_StatusTest_FiniteValue.C \
	$(srcdir)/NOX_StatusTest_NormUpdate.C \
	$(srcdir)/NOX_StatusTest_Combo.C \
	$(srcdir)/NOX_StatusTest_MaxIters.C \
        $(srcdir)/NOX_StatusTest_NormWRMS.C \
        $(srcdir)/NOX_StatusTest_Stagnation.C \
        $(srcdir)/NOX_StatusTest_Divergence.C 

UTILS_H = \
	$(srcdir)/NOX_Utils.H \
        $(srcdir)/NOX_SharedObjectTemplate.H \
	$(srcdir)/NOX_Random.H

UTILS = \
	$(srcdir)/NOX_Utils.C \
	$(srcdir)/NOX_Random.C

if BUILD_NOX_MULTIVECS

MULTIVECS_H = \
	$(srcdir)/NOX_Abstract_MultiVector.H \
	$(srcdir)/NOX_MultiVector.H

MULTIVECS = \
	$(srcdir)/NOX_MultiVector.C

else

MULTIVECS_H =
MULTIVECS =

endif


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

libnox_a_LIBADD = $(CXXREPOS) $(XTRALDADD)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Files to be deleted by 'make clean'
# ------------------------------------------------------------------------

MOSTLYCLEANFILES = $(CXXREPOS)





