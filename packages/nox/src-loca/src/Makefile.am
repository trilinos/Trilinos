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


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside NOX. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

LOCA_OLD_H = \
	$(srcdir)/LOCA_SingularJacobianSolve_Generic.H \
	$(srcdir)/LOCA_SingularJacobianSolve_Manager.H \
	$(srcdir)/LOCA_SingularJacobianSolve_Default.H \
	$(srcdir)/LOCA_SingularJacobianSolve_Nic.H \
	$(srcdir)/LOCA_SingularJacobianSolve_NicDay.H \
	$(srcdir)/LOCA_SingularJacobianSolve_ItRef.H \
	$(srcdir)/LOCA_Continuation_StatusTest_ParameterUpdateNorm.H \
	$(srcdir)/LOCA_Continuation_StatusTest_ParameterResidualNorm.H \
	$(srcdir)/LOCA_Bifurcation_TPBord_StatusTest_ParameterUpdateNorm.H \
	$(srcdir)/LOCA_Bifurcation_TPBord_StatusTest_NullVectorNormWRMS.H \
	$(srcdir)/LOCA_Bifurcation_PitchforkBord_ParameterUpdateNorm.H \
	$(srcdir)/LOCA_Bifurcation_PitchforkBord_SlackUpdateNorm.H \
	$(srcdir)/LOCA_Bifurcation_PitchforkBord_NullVectorNormWRMS.H

LOCA_OLD = \
	$(srcdir)/LOCA_SingularJacobianSolve_Manager.C \
	$(srcdir)/LOCA_SingularJacobianSolve_Default.C \
	$(srcdir)/LOCA_SingularJacobianSolve_Nic.C \
	$(srcdir)/LOCA_SingularJacobianSolve_NicDay.C \
	$(srcdir)/LOCA_SingularJacobianSolve_ItRef.C \
	$(srcdir)/LOCA_Continuation_StatusTest_ParameterUpdateNorm.C \
	$(srcdir)/LOCA_Continuation_StatusTest_ParameterResidualNorm.C \
	$(srcdir)/LOCA_Bifurcation_TPBord_StatusTest_ParameterUpdateNorm.C \
	$(srcdir)/LOCA_Bifurcation_TPBord_StatusTest_NullVectorNormWRMS.C \
	$(srcdir)/LOCA_Bifurcation_PitchforkBord_ParameterUpdateNorm.C \
	$(srcdir)/LOCA_Bifurcation_PitchforkBord_SlackUpdateNorm.C \
	$(srcdir)/LOCA_Bifurcation_PitchforkBord_NullVectorNormWRMS.C 

LOCA_ABSTRACT_H = \
        $(srcdir)/LOCA_Abstract_Group.H \
	$(srcdir)/LOCA_Abstract_TransposeSolveGroup.H \
	$(srcdir)/LOCA_Abstract_Iterator.H \
	$(srcdir)/LOCA_Abstract_Factory.H

LOCA_ABSTRACT = \
        $(srcdir)/LOCA_Abstract_Group.C \
	$(srcdir)/LOCA_Abstract_Iterator.C \
	$(srcdir)/LOCA_Abstract_Factory.C

LOCA_EXTENDED_H = \
	$(srcdir)/LOCA_Extended_Vector.H \
	$(srcdir)/LOCA_Extended_MultiVector.H \
	$(srcdir)/LOCA_Extended_MultiAbstractGroup.H

LOCA_EXTENDED = \
	$(srcdir)/LOCA_Extended_Vector.C \
	$(srcdir)/LOCA_Extended_MultiVector.C \
	$(srcdir)/LOCA_Extended_MultiAbstractGroup.C

LOCA_BORDERED_SOLVER_H = \
	$(srcdir)/LOCA_BorderedSolver_AbstractStrategy.H \
	$(srcdir)/LOCA_BorderedSolver_LowerTriangularBlockElimination.H \
	$(srcdir)/LOCA_BorderedSolver_UpperTriangularBlockElimination.H \
	$(srcdir)/LOCA_BorderedSolver_Bordering.H \
	$(srcdir)/LOCA_BorderedSolver_HouseholderQR.H \
	$(srcdir)/LOCA_BorderedSystem_AbstractGroup.H \
	$(srcdir)/LOCA_BorderedSolver_Nested.H \
	$(srcdir)/LOCA_BorderedSolver_Factory.H

LOCA_BORDERED_SOLVER = \
	$(srcdir)/LOCA_BorderedSolver_LowerTriangularBlockElimination.C \
	$(srcdir)/LOCA_BorderedSolver_UpperTriangularBlockElimination.C \
	$(srcdir)/LOCA_BorderedSolver_AbstractStrategy.C \
	$(srcdir)/LOCA_BorderedSolver_Bordering.C \
	$(srcdir)/LOCA_BorderedSolver_HouseholderQR.C \
	$(srcdir)/LOCA_BorderedSolver_Nested.C \
	$(srcdir)/LOCA_BorderedSolver_Factory.C

LOCA_MULTICONTINUATION_H = \
	$(srcdir)/LOCA_MultiContinuation_ExtendedVector.H \
	$(srcdir)/LOCA_MultiContinuation_ExtendedMultiVector.H \
	$(srcdir)/LOCA_MultiContinuation_AbstractGroup.H \
	$(srcdir)/LOCA_MultiContinuation_FiniteDifferenceGroup.H \
	$(srcdir)/LOCA_MultiContinuation_ConstraintInterface.H \
	$(srcdir)/LOCA_MultiContinuation_ConstraintInterfaceMVDX.H \
	$(srcdir)/LOCA_MultiContinuation_MultiVecConstraint.H \
	$(srcdir)/LOCA_MultiContinuation_CompositeConstraint.H \
	$(srcdir)/LOCA_MultiContinuation_CompositeConstraintMVDX.H \
	$(srcdir)/LOCA_MultiContinuation_ConstrainedGroup.H \
	$(srcdir)/LOCA_MultiContinuation_AbstractStrategy.H \
	$(srcdir)/LOCA_MultiContinuation_ExtendedGroup.H \
	$(srcdir)/LOCA_MultiContinuation_NaturalGroup.H \
	$(srcdir)/LOCA_MultiContinuation_NaturalConstraint.H \
	$(srcdir)/LOCA_MultiContinuation_ArcLengthGroup.H \
	$(srcdir)/LOCA_MultiContinuation_ArcLengthConstraint.H \
	$(srcdir)/LOCA_MultiContinuation_Factory.H 

LOCA_MULTICONTINUATION = \
	$(srcdir)/LOCA_MultiContinuation_ExtendedVector.C \
	$(srcdir)/LOCA_MultiContinuation_ExtendedMultiVector.C \
	$(srcdir)/LOCA_MultiContinuation_AbstractGroup.C \
	$(srcdir)/LOCA_MultiContinuation_FiniteDifferenceGroup.C \
	$(srcdir)/LOCA_MultiContinuation_ConstraintInterfaceMVDX.C \
	$(srcdir)/LOCA_MultiContinuation_MultiVecConstraint.C \
	$(srcdir)/LOCA_MultiContinuation_CompositeConstraint.C \
	$(srcdir)/LOCA_MultiContinuation_CompositeConstraintMVDX.C \
	$(srcdir)/LOCA_MultiContinuation_ConstrainedGroup.C \
	$(srcdir)/LOCA_MultiContinuation_ExtendedGroup.C \
	$(srcdir)/LOCA_MultiContinuation_NaturalGroup.C \
	$(srcdir)/LOCA_MultiContinuation_NaturalConstraint.C \
	$(srcdir)/LOCA_MultiContinuation_ArcLengthGroup.C \
	$(srcdir)/LOCA_MultiContinuation_ArcLengthConstraint.C \
	$(srcdir)/LOCA_MultiContinuation_Factory.C

LOCA_SINGULARSOLVE_H = 

LOCA_SINGULARSOLVE = 

LOCA_BIFURCATION_H = \
	$(srcdir)/LOCA_TimeDependent_AbstractGroup.H \
	$(srcdir)/LOCA_TurningPoint_MooreSpence_ExtendedVector.H \
	$(srcdir)/LOCA_TurningPoint_MooreSpence_ExtendedMultiVector.H \
	$(srcdir)/LOCA_TurningPoint_MooreSpence_AbstractGroup.H \
	$(srcdir)/LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H \
	$(srcdir)/LOCA_TurningPoint_MooreSpence_ExtendedGroup.H \
	$(srcdir)/LOCA_TurningPoint_MooreSpence_SolverStrategy.H \
	$(srcdir)/LOCA_TurningPoint_MooreSpence_SalingerBordering.H \
	$(srcdir)/LOCA_TurningPoint_MooreSpence_PhippsBordering.H \
	$(srcdir)/LOCA_TurningPoint_MooreSpence_SolverFactory.H \
	$(srcdir)/LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H \
	$(srcdir)/LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.H\
	$(srcdir)/LOCA_TurningPoint_MinimallyAugmented_Constraint.H \
	$(srcdir)/LOCA_TurningPoint_MinimallyAugmented_ModifiedConstraint.H \
	$(srcdir)/LOCA_TurningPoint_MinimallyAugmented_ExtendedGroup.H \
	$(srcdir)/LOCA_Pitchfork_MooreSpence_ExtendedVector.H \
	$(srcdir)/LOCA_Pitchfork_MooreSpence_ExtendedMultiVector.H \
	$(srcdir)/LOCA_Pitchfork_MooreSpence_AbstractGroup.H \
	$(srcdir)/LOCA_Pitchfork_MooreSpence_ExtendedGroup.H \
	$(srcdir)/LOCA_Pitchfork_MooreSpence_SolverStrategy.H \
	$(srcdir)/LOCA_Pitchfork_MooreSpence_SalingerBordering.H \
	$(srcdir)/LOCA_Pitchfork_MooreSpence_PhippsBordering.H \
	$(srcdir)/LOCA_Pitchfork_MooreSpence_SolverFactory.H \
	$(srcdir)/LOCA_Hopf_MooreSpence_ExtendedVector.H \
	$(srcdir)/LOCA_Hopf_MooreSpence_ExtendedMultiVector.H \
	$(srcdir)/LOCA_Hopf_MooreSpence_AbstractGroup.H \
	$(srcdir)/LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.H \
	$(srcdir)/LOCA_Hopf_MooreSpence_ExtendedGroup.H \
	$(srcdir)/LOCA_Hopf_MooreSpence_SolverStrategy.H \
	$(srcdir)/LOCA_Hopf_MooreSpence_SalingerBordering.H \
	$(srcdir)/LOCA_Hopf_MooreSpence_SolverFactory.H \
	$(srcdir)/LOCA_Bifurcation_Factory.H

LOCA_BIFURCATION = \
	$(srcdir)/LOCA_TurningPoint_MooreSpence_ExtendedVector.C \
	$(srcdir)/LOCA_TurningPoint_MooreSpence_ExtendedMultiVector.C \
	$(srcdir)/LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.C \
	$(srcdir)/LOCA_TurningPoint_MooreSpence_ExtendedGroup.C \
	$(srcdir)/LOCA_TurningPoint_MooreSpence_SalingerBordering.C \
	$(srcdir)/LOCA_TurningPoint_MooreSpence_PhippsBordering.C \
	$(srcdir)/LOCA_TurningPoint_MooreSpence_SolverFactory.C \
	$(srcdir)/LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.C\
	$(srcdir)/LOCA_TurningPoint_MinimallyAugmented_Constraint.C \
	$(srcdir)/LOCA_TurningPoint_MinimallyAugmented_ModifiedConstraint.C \
	$(srcdir)/LOCA_TurningPoint_MinimallyAugmented_ExtendedGroup.C \
	$(srcdir)/LOCA_Pitchfork_MooreSpence_ExtendedVector.C \
	$(srcdir)/LOCA_Pitchfork_MooreSpence_ExtendedMultiVector.C \
	$(srcdir)/LOCA_Pitchfork_MooreSpence_ExtendedGroup.C \
	$(srcdir)/LOCA_Pitchfork_MooreSpence_SalingerBordering.C \
	$(srcdir)/LOCA_Pitchfork_MooreSpence_PhippsBordering.C \
	$(srcdir)/LOCA_Pitchfork_MooreSpence_SolverFactory.C \
	$(srcdir)/LOCA_Hopf_MooreSpence_ExtendedVector.C \
	$(srcdir)/LOCA_Hopf_MooreSpence_ExtendedMultiVector.C \
	$(srcdir)/LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.C \
	$(srcdir)/LOCA_Hopf_MooreSpence_ExtendedGroup.C \
	$(srcdir)/LOCA_Hopf_MooreSpence_SalingerBordering.C \
	$(srcdir)/LOCA_Hopf_MooreSpence_SolverFactory.C \
	$(srcdir)/LOCA_Bifurcation_Factory.C

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
        $(srcdir)/LOCA_ErrorCheck.H

LOCA_UTILS = \
        $(srcdir)/LOCA_ErrorCheck.C

LOCA_PARAMETER_H = \
        $(srcdir)/LOCA_Parameter_Vector.H \
	$(srcdir)/LOCA_Parameter_Library.H \
	$(srcdir)/LOCA_Parameter_LibraryT.H \
	$(srcdir)/LOCA_Parameter_Entry.H \
	$(srcdir)/LOCA_Parameter_SublistParser.H

LOCA_PARAMETER = \
        $(srcdir)/LOCA_Parameter_Vector.C \
	$(srcdir)/LOCA_Parameter_Library.C \
	$(srcdir)/LOCA_Parameter_SublistParser.C

LOCA_PREDICTOR_H = \
	$(srcdir)/LOCA_MultiPredictor_AbstractStrategy.H \
	$(srcdir)/LOCA_MultiPredictor_Constant.H \
	$(srcdir)/LOCA_MultiPredictor_Tangent.H \
	$(srcdir)/LOCA_MultiPredictor_Secant.H \
	$(srcdir)/LOCA_MultiPredictor_Random.H \
	$(srcdir)/LOCA_MultiPredictor_Restart.H \
	$(srcdir)/LOCA_MultiPredictor_Factory.H

LOCA_PREDICTOR = \
	$(srcdir)/LOCA_MultiPredictor_AbstractStrategy.C \
	$(srcdir)/LOCA_MultiPredictor_Constant.C \
	$(srcdir)/LOCA_MultiPredictor_Tangent.C \
	$(srcdir)/LOCA_MultiPredictor_Secant.C \
	$(srcdir)/LOCA_MultiPredictor_Random.C \
	$(srcdir)/LOCA_MultiPredictor_Restart.C \
	$(srcdir)/LOCA_MultiPredictor_Factory.C

LOCA_STEPSIZE_H = \
	$(srcdir)/LOCA_StepSize_AbstractStrategy.H \
	$(srcdir)/LOCA_StepSize_Constant.H \
	$(srcdir)/LOCA_StepSize_Adaptive.H \
	$(srcdir)/LOCA_StepSize_Factory.H

LOCA_STEPSIZE = \
	$(srcdir)/LOCA_StepSize_Constant.C \
	$(srcdir)/LOCA_StepSize_Adaptive.C \
	$(srcdir)/LOCA_StepSize_Factory.C

LOCA_STEPPER_H = \
	$(srcdir)/LOCA_Stepper.H

LOCA_STEPPER = \
	$(srcdir)/LOCA_Stepper.C

LOCA_STATUSTEST_H = \
	$(srcdir)/LOCA_Solver_Wrapper.H \
	$(srcdir)/LOCA_StatusTest_Wrapper.H

LOCA_STATUSTEST = \
	$(srcdir)/LOCA_Solver_Wrapper.C \
	$(srcdir)/LOCA_StatusTest_Wrapper.C

LOCA_EIGENSOLVER_H = \
	$(srcdir)/LOCA_Eigensolver_AbstractStrategy.H \
	$(srcdir)/LOCA_Eigensolver_DefaultStrategy.H \
	$(srcdir)/LOCA_Eigensolver_Factory.H \
	$(srcdir)/LOCA_EigenvalueSort_Strategies.H \
	$(srcdir)/LOCA_EigenvalueSort_Factory.H \
	$(srcdir)/LOCA_SaveEigenData_AbstractStrategy.H \
	$(srcdir)/LOCA_SaveEigenData_DefaultStrategy.H \
	$(srcdir)/LOCA_SaveEigenData_Factory.H

LOCA_EIGENSOLVER = \
	$(srcdir)/LOCA_Eigensolver_DefaultStrategy.C \
	$(srcdir)/LOCA_Eigensolver_Factory.C \
	$(srcdir)/LOCA_EigenvalueSort_Strategies.C \
	$(srcdir)/LOCA_EigenvalueSort_Factory.C \
	$(srcdir)/LOCA_SaveEigenData_DefaultStrategy.C \
	$(srcdir)/LOCA_SaveEigenData_Factory.C

LOCA_FACTORY_H = \
	$(srcdir)/LOCA_Factory.H \
	$(srcdir)/LOCA_GlobalData.H 

LOCA_FACTORY = \
	$(srcdir)/LOCA_Factory.C \
	$(srcdir)/LOCA_GlobalData.C

EXTRA_LOCA_ANASAZI_H = \
	$(srcdir)/Anasazi_LOCA_MultiVecTraits.H \
	$(srcdir)/Anasazi_LOCA_OperatorTraits.H \
	$(srcdir)/Anasazi_LOCA_Sort.H \
	$(srcdir)/LOCA_Eigensolver_AnasaziStrategy.H

EXTRA_LOCA_ANASAZI = \
	$(srcdir)/Anasazi_LOCA_Sort.C \
	$(srcdir)/LOCA_Eigensolver_AnasaziStrategy.C

if BUILD_LOCAANASAZI

LOCA_ANASAZI_H = \
	$(srcdir)/LOCA_AnasaziOperator_AbstractStrategy.H \
	$(srcdir)/LOCA_AnasaziOperator_JacobianInverse.H \
	$(srcdir)/LOCA_AnasaziOperator_ShiftInvert.H \
	$(srcdir)/LOCA_AnasaziOperator_Cayley.H \
	$(srcdir)/LOCA_AnasaziOperator_Factory.H \
	$(EXTRA_LOCA_ANASAZI_H)

LOCA_ANASAZI = \
	$(srcdir)/LOCA_AnasaziOperator_JacobianInverse.C \
	$(srcdir)/LOCA_AnasaziOperator_ShiftInvert.C \
	$(srcdir)/LOCA_AnasaziOperator_Cayley.C \
	$(srcdir)/LOCA_AnasaziOperator_Factory.C \
	$(EXTRA_LOCA_ANASAZI)

else

LOCA_ANASAZI_H = \
	$(srcdir)/LOCA_AnasaziOperator_AbstractStrategy.H \
	$(srcdir)/LOCA_AnasaziOperator_JacobianInverse.H \
	$(srcdir)/LOCA_AnasaziOperator_ShiftInvert.H \
	$(srcdir)/LOCA_AnasaziOperator_Cayley.H \
	$(srcdir)/LOCA_AnasaziOperator_Factory.H

LOCA_ANASAZI = \
	$(srcdir)/LOCA_AnasaziOperator_JacobianInverse.C \
	$(srcdir)/LOCA_AnasaziOperator_ShiftInvert.C \
	$(srcdir)/LOCA_AnasaziOperator_Cayley.C \
	$(srcdir)/LOCA_AnasaziOperator_Factory.C

endif

LOCA_H = \
        $(srcdir)/LOCA.H \
	$(LOCA_EXTENDED_H) \
	$(LOCA_BORDERED_SOLVER_H) \
	$(LOCA_MULTICONTINUATION_H) \
	$(LOCA_ANASAZI_H) \
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
	$(LOCA_STATUSTEST_H) \
	$(LOCA_EIGENSOLVER_H) \
	$(LOCA_FACTORY_H)

LOCA = \
	$(LOCA_EXTENDED) \
	$(LOCA_BORDERED_SOLVER) \
	$(LOCA_MULTICONTINUATION) \
	$(LOCA_ANASAZI) \
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
	$(LOCA_STATUSTEST) \
	$(LOCA_EIGENSOLVER) \
	$(LOCA_FACTORY)

# ------------------------------------------------------------------------
# Includes for required and optional library dependencies
# ------------------------------------------------------------------------
include $(top_builddir)/Makefile.export.nox

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(LOCA_INCLUDES))
else
EXPORT_INCLUDES = $(LOCA_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# ------------------------------------------------------------------------
# Loca library specifications
# ------------------------------------------------------------------------

libloca_a_SOURCES = $(LOCA_H) $(LOCA)

EXTRA_libloca_a_SOURCES = $(EXTRA_LOCA_ANASAZI) $(EXTRA_LOCA_ANASAZI_H) 

libloca_a_CXXFLAGS = \
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/src-lapack \
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
