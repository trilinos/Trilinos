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
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside NOX. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

LOCA_ABSTRACT_H = \
        $(srcdir)/LOCA_Abstract_Group.H \
	$(srcdir)/LOCA_Abstract_Iterator.H

LOCA_ABSTRACT = \
        $(srcdir)/LOCA_Abstract_Group.C \
	$(srcdir)/LOCA_Abstract_Iterator.C

LOCA_EXTENDED_H = \
	$(srcdir)/LOCA_Extended_Vector.H \
	$(srcdir)/LOCA_Extended_MultiVector.H \
	$(srcdir)/LOCA_Extended_AbstractGroup.H

LOCA_EXTENDED = \
	$(srcdir)/LOCA_Extended_Vector.C \
	$(srcdir)/LOCA_Extended_MultiVector.C \
	$(srcdir)/LOCA_Extended_AbstractGroup.C

LOCA_BORDERED_SYSTEM_H = \
	$(srcdir)/LOCA_BorderedSystem_Generic.H \
	$(srcdir)/LOCA_BorderedSystem_Manager.H \
	$(srcdir)/LOCA_BorderedSystem_Bordering.H

LOCA_BORDERED_SYSTEM = \
	$(srcdir)/LOCA_BorderedSystem_Manager.C \
	$(srcdir)/LOCA_BorderedSystem_Bordering.C

LOCA_CONTINUATION_H = \
	$(srcdir)/LOCA_Continuation_ExtendedVector.H \
	$(srcdir)/LOCA_Continuation_AbstractGroup.H \
	$(srcdir)/LOCA_Continuation_FiniteDifferenceGroup.H \
	$(srcdir)/LOCA_Continuation_AnasaziGroup.H \
	$(srcdir)/LOCA_Continuation_ExtendedGroup.H \
	$(srcdir)/LOCA_Continuation_Manager.H \
	$(srcdir)/LOCA_Continuation_NaturalGroup.H \
	$(srcdir)/LOCA_Continuation_ArcLengthGroup.H \
	$(srcdir)/LOCA_Continuation_HouseholderGroup.H

LOCA_CONTINUATION = \
	$(srcdir)/LOCA_Continuation_ExtendedVector.C \
	$(srcdir)/LOCA_Continuation_AbstractGroup.C \
	$(srcdir)/LOCA_Continuation_AnasaziGroup.C \
	$(srcdir)/LOCA_Continuation_FiniteDifferenceGroup.C \
	$(srcdir)/LOCA_Continuation_ExtendedGroup.C \
	$(srcdir)/LOCA_Continuation_Manager.C \
	$(srcdir)/LOCA_Continuation_NaturalGroup.C \
	$(srcdir)/LOCA_Continuation_ArcLengthGroup.C \
	$(srcdir)/LOCA_Continuation_HouseholderGroup.C

LOCA_MULTICONTINUATION_H = \
	$(srcdir)/LOCA_MultiContinuation_ExtendedVector.H \
	$(srcdir)/LOCA_MultiContinuation_ExtendedMultiVector.H \
	$(srcdir)/LOCA_MultiContinuation_AbstractGroup.H \
	$(srcdir)/LOCA_MultiContinuation_FiniteDifferenceGroup.H \
	$(srcdir)/LOCA_MultiContinuation_ExtendedGroup.H \
	$(srcdir)/LOCA_MultiContinuation_ArcLengthGroup.H \
	$(srcdir)/LOCA_MultiContinuation_NaturalGroup.H \
	$(srcdir)/LOCA_MultiContinuation_Manager.H

LOCA_MULTICONTINUATION = \
	$(srcdir)/LOCA_MultiContinuation_ExtendedVector.C \
	$(srcdir)/LOCA_MultiContinuation_ExtendedMultiVector.C \
	$(srcdir)/LOCA_MultiContinuation_FiniteDifferenceGroup.C \
	$(srcdir)/LOCA_MultiContinuation_ExtendedGroup.C \
	$(srcdir)/LOCA_MultiContinuation_ArcLengthGroup.C \
	$(srcdir)/LOCA_MultiContinuation_NaturalGroup.C \
	$(srcdir)/LOCA_MultiContinuation_Manager.C

LOCA_SINGULARSOLVE_H = \
	$(srcdir)/LOCA_SingularJacobianSolve_Generic.H \
	$(srcdir)/LOCA_SingularJacobianSolve_Manager.H \
	$(srcdir)/LOCA_SingularJacobianSolve_Default.H \
	$(srcdir)/LOCA_SingularJacobianSolve_Nic.H \
	$(srcdir)/LOCA_SingularJacobianSolve_NicDay.H \
	$(srcdir)/LOCA_SingularJacobianSolve_ItRef.H

LOCA_SINGULARSOLVE = \
	$(srcdir)/LOCA_SingularJacobianSolve_Manager.C \
	$(srcdir)/LOCA_SingularJacobianSolve_Default.C \
	$(srcdir)/LOCA_SingularJacobianSolve_Nic.C \
	$(srcdir)/LOCA_SingularJacobianSolve_NicDay.C \
	$(srcdir)/LOCA_SingularJacobianSolve_ItRef.C

LOCA_BIFURCATION_H = \
        $(srcdir)/LOCA_Bifurcation_TPBord_ExtendedVector.H \
	$(srcdir)/LOCA_Bifurcation_TPBord_AbstractGroup.H \
	$(srcdir)/LOCA_Bifurcation_TPBord_FiniteDifferenceGroup.H \
	$(srcdir)/LOCA_Bifurcation_TPBord_SingularSolveGroup.H \
        $(srcdir)/LOCA_Bifurcation_TPBord_ExtendedGroup.H \
	$(srcdir)/LOCA_Bifurcation_TPBord_ModifiedBorderingGroup.H \
	$(srcdir)/LOCA_Bifurcation_TPBord_NicDayModifiedBorderingGroup.H \
	$(srcdir)/LOCA_Bifurcation_PitchforkBord_ExtendedVector.H \
	$(srcdir)/LOCA_Bifurcation_PitchforkBord_ExtendedGroup.H \
	$(srcdir)/LOCA_TimeDependent_AbstractGroup.H \
	$(srcdir)/LOCA_Bifurcation_HopfBord_ExtendedVector.H \
	$(srcdir)/LOCA_Bifurcation_HopfBord_AbstractGroup.H \
	$(srcdir)/LOCA_Bifurcation_HopfBord_FiniteDifferenceGroup.H \
	$(srcdir)/LOCA_Bifurcation_HopfBord_ExtendedGroup.H \
	$(srcdir)/LOCA_Bifurcation_Manager.H

LOCA_BIFURCATION = \
        $(srcdir)/LOCA_Bifurcation_TPBord_ExtendedVector.C \
	$(srcdir)/LOCA_Bifurcation_TPBord_AbstractGroup.C \
	$(srcdir)/LOCA_Bifurcation_TPBord_FiniteDifferenceGroup.C \
	$(srcdir)/LOCA_Bifurcation_TPBord_SingularSolveGroup.C \
	$(srcdir)/LOCA_Bifurcation_TPBord_ExtendedGroup.C \
	$(srcdir)/LOCA_Bifurcation_TPBord_ModifiedBorderingGroup.C \
	$(srcdir)/LOCA_Bifurcation_TPBord_NicDayModifiedBorderingGroup.C \
	$(srcdir)/LOCA_Bifurcation_PitchforkBord_ExtendedVector.C \
	$(srcdir)/LOCA_Bifurcation_PitchforkBord_ExtendedGroup.C \
	$(srcdir)/LOCA_Bifurcation_HopfBord_ExtendedVector.C \
	$(srcdir)/LOCA_Bifurcation_HopfBord_AbstractGroup.C \
	$(srcdir)/LOCA_Bifurcation_HopfBord_FiniteDifferenceGroup.C \
	$(srcdir)/LOCA_Bifurcation_HopfBord_ExtendedGroup.C \
	$(srcdir)/LOCA_Bifurcation_Manager.C

LOCA_HOMOTOPY_H = \
	$(srcdir)/LOCA_Homotopy_AbstractGroup.H \
        $(srcdir)/LOCA_Homotopy_Group.H

LOCA_HOMOTOPY = \
	$(srcdir)/LOCA_Homotopy_Group.C

LOCA_DERIVUTILS_H = \
        $(srcdir)/LOCA_DerivUtils.H

LOCA_DERIVUTILS = \
        $(srcdir)/LOCA_DerivUtils.C

LOCA_UTILS_H = \
        $(srcdir)/LOCA_Utils.H \
        $(srcdir)/LOCA_ErrorCheck.H

LOCA_UTILS = \
        $(srcdir)/LOCA_Utils.C \
        $(srcdir)/LOCA_ErrorCheck.C

LOCA_PARAMETER_H = \
        $(srcdir)/LOCA_Parameter_Vector.H \
	$(srcdir)/LOCA_Parameter_Library.H \
	$(srcdir)/LOCA_Parameter_LibraryT.H \
	$(srcdir)/LOCA_Parameter_Entry.H

LOCA_PARAMETER = \
        $(srcdir)/LOCA_Parameter_Vector.C \
	$(srcdir)/LOCA_Parameter_Library.C 

LOCA_PREDICTOR_H = \
	$(srcdir)/LOCA_Predictor_Generic.H \
	$(srcdir)/LOCA_Predictor_Manager.H \
	$(srcdir)/LOCA_Predictor_Constant.H \
	$(srcdir)/LOCA_Predictor_Tangent.H \
	$(srcdir)/LOCA_Predictor_Secant.H \
	$(srcdir)/LOCA_Predictor_Random.H

LOCA_PREDICTOR = \
	$(srcdir)/LOCA_Predictor_Generic.C \
	$(srcdir)/LOCA_Predictor_Manager.C \
	$(srcdir)/LOCA_Predictor_Constant.C \
	$(srcdir)/LOCA_Predictor_Tangent.C \
	$(srcdir)/LOCA_Predictor_Secant.C \
	$(srcdir)/LOCA_Predictor_Random.C

LOCA_STEPSIZE_H = \
	$(srcdir)/LOCA_StepSize_Generic.H \
	$(srcdir)/LOCA_StepSize_Manager.H \
	$(srcdir)/LOCA_StepSize_Constant.H \
	$(srcdir)/LOCA_StepSize_Adaptive.H 

LOCA_STEPSIZE = \
	$(srcdir)/LOCA_StepSize_Manager.C \
	$(srcdir)/LOCA_StepSize_Constant.C \
	$(srcdir)/LOCA_StepSize_Adaptive.C 

LOCA_STEPPER_H = \
	$(srcdir)/LOCA_Stepper.H \
	$(srcdir)/LOCA_NewStepper.H

LOCA_STEPPER = \
	$(srcdir)/LOCA_Stepper.C \
	$(srcdir)/LOCA_NewStepper.C

LOCA_STATUSTEST_H = \
	$(srcdir)/LOCA_Solver_Wrapper.H \
	$(srcdir)/LOCA_StatusTest_Wrapper.H \
	$(srcdir)/LOCA_Continuation_StatusTest_ParameterUpdateNorm.H \
	$(srcdir)/LOCA_Continuation_StatusTest_ParameterResidualNorm.H \
	$(srcdir)/LOCA_Bifurcation_TPBord_StatusTest_ParameterUpdateNorm.H \
	$(srcdir)/LOCA_Bifurcation_TPBord_StatusTest_NullVectorNormWRMS.H \
	$(srcdir)/LOCA_Bifurcation_PitchforkBord_StatusTest_ParameterUpdateNorm.H \
	$(srcdir)/LOCA_Bifurcation_PitchforkBord_StatusTest_SlackUpdateNorm.H \
	$(srcdir)/LOCA_Bifurcation_PitchforkBord_StatusTest_NullVectorNormWRMS.H

LOCA_STATUSTEST = \
	$(srcdir)/LOCA_Solver_Wrapper.C \
	$(srcdir)/LOCA_StatusTest_Wrapper.C \
	$(srcdir)/LOCA_Continuation_StatusTest_ParameterUpdateNorm.C \
	$(srcdir)/LOCA_Continuation_StatusTest_ParameterResidualNorm.C \
	$(srcdir)/LOCA_Bifurcation_TPBord_StatusTest_ParameterUpdateNorm.C \
	$(srcdir)/LOCA_Bifurcation_TPBord_StatusTest_NullVectorNormWRMS.C \
	$(srcdir)/LOCA_Bifurcation_PitchforkBord_StatusTest_ParameterUpdateNorm.C \
	$(srcdir)/LOCA_Bifurcation_PitchforkBord_StatusTest_SlackUpdateNorm.C \
	$(srcdir)/LOCA_Bifurcation_PitchforkBord_StatusTest_NullVectorNormWRMS.C

LOCA_H = \
        $(srcdir)/LOCA.H \
	$(LOCA_EXTENDED_H) \
	$(LOCA_BORDERED_SYSTEM_H) \
	$(LOCA_CONTINUATION_H) \
	$(LOCA_MULTICONTINUATION_H) \
	$(LOCA_SINGULARSOLVE_H) \
	$(LOCA_BIFURCATION_H) \
	$(LOCA_HOMOTOPY_H) \
	$(LOCA_DERIVUTILS_H) \
	$(LOCA_UTILS_H) \
	$(LOCA_PARAMETER_H) \
	$(LOCA_PREDICTOR_H) \
	$(LOCA_STEPSIZE_H) \
	$(LOCA_STEPPER_H) \
	$(LOCA_ABSTRACT_H) \
	$(LOCA_STATUSTEST_H)

LOCA = \
	$(LOCA_EXTENDED) \
	$(LOCA_BORDERED_SYSTEM) \
	$(LOCA_CONTINUATION) \
	$(LOCA_MULTICONTINUATION) \
	$(LOCA_SINGULARSOLVE) \
	$(LOCA_BIFURCATION) \
	$(LOCA_HOMOTOPY) \
	$(LOCA_DERIVUTILS) \
	$(LOCA_UTILS) \
	$(LOCA_PARAMETER) \
	$(LOCA_PREDICTOR) \
	$(LOCA_STEPSIZE) \
	$(LOCA_STEPPER) \
	$(LOCA_ABSTRACT) \
	$(LOCA_STATUSTEST)

# ------------------------------------------------------------------------
# Loca library specifications
# ------------------------------------------------------------------------

libloca_a_SOURCES = $(LOCA_H) $(LOCA)

libloca_a_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src-lapack \
                      $(TRILINOS_CXXFLAGS) 

# ------------------------------------------------------------------------
# When to turn it on
# ------------------------------------------------------------------------

if BUILD_LIBLOCA

# List of all libraries to install in $(libexecdir)

lib_LIBRARIES = libloca.a

# List of all headers to install in $(includedir)

include_HEADERS = $(LOCA_H)

endif

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libloca_a_AR = $(ALTERNATE_AR)

else

libloca_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libloca_a_LIBADD = $(CXXREPOS) $(XTRALDADD)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Files to be deleted by 'make clean'
# ------------------------------------------------------------------------

MOSTLYCLEANFILES = $(CXXREPOS)

