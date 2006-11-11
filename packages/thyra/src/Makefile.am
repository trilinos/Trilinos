# @HEADER
# ************************************************************************
# 
#  Thyra: Interfaces and Support Code for the Interoperability of
#                Abstract Numerical Algorithms
#
#                Copyright (2004) Sandia Corporation
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
	@echo "Trilinos package thyra subdirectory src built successfully."
	@echo ""

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
# There is, at present, no distinction between most of the categories, 
# other than those which start with EXTRA_.  All categories which do not
# start with EXTRA_ are included in all builds.  Categories which start
# with EXTRA_ are included in some builds (typically controlled by 
# --with-xyz or --without-xyz)
#

INTERFACES_OPERATOR_VECTOR_ANA_FUNDAMENTAL_H = \
  $(srcdir)/interfaces/operator_vector/ana/fundamental/Thyra_LinearOpBaseDecl.hpp\
  $(srcdir)/interfaces/operator_vector/ana/fundamental/Thyra_LinearOpBase.hpp\
  $(srcdir)/interfaces/operator_vector/ana/fundamental/Thyra_MultiVectorBaseDecl.hpp\
  $(srcdir)/interfaces/operator_vector/ana/fundamental/Thyra_MultiVectorBase.hpp\
  $(srcdir)/interfaces/operator_vector/ana/fundamental/Thyra_OperatorVectorTypes.hpp\
  $(srcdir)/interfaces/operator_vector/ana/fundamental/Thyra_VectorBaseDecl.hpp\
  $(srcdir)/interfaces/operator_vector/ana/fundamental/Thyra_VectorBase.hpp\
  $(srcdir)/interfaces/operator_vector/ana/fundamental/Thyra_VectorSpaceBaseDecl.hpp\
  $(srcdir)/interfaces/operator_vector/ana/fundamental/Thyra_VectorSpaceBase.hpp\
  $(srcdir)/interfaces/operator_vector/ana/fundamental/Thyra_VectorSpaceFactoryBaseDecl.hpp\
  $(srcdir)/interfaces/operator_vector/ana/fundamental/Thyra_VectorSpaceFactoryBase.hpp

INTERFACES_OPERATOR_VECTOR_ANA_EXTENDED_H = \
  $(srcdir)/interfaces/operator_vector/ana/extended/Thyra_AddedLinearOpBase.hpp\
  $(srcdir)/interfaces/operator_vector/ana/extended/Thyra_BlockedLinearOpBase.hpp\
  $(srcdir)/interfaces/operator_vector/ana/extended/Thyra_DiagonalLinearOpBase.hpp\
  $(srcdir)/interfaces/operator_vector/ana/extended/Thyra_IdentityLinearOpBase.hpp\
  $(srcdir)/interfaces/operator_vector/ana/extended/Thyra_MultipliedLinearOpBase.hpp\
  $(srcdir)/interfaces/operator_vector/ana/extended/Thyra_MultiVectorRandomizerBase.hpp\
  $(srcdir)/interfaces/operator_vector/ana/extended/Thyra_PhysicallyBlockedLinearOpBase.hpp\
  $(srcdir)/interfaces/operator_vector/ana/extended/Thyra_PolynomialVectorTraits.hpp\
  $(srcdir)/interfaces/operator_vector/ana/extended/Thyra_ProductMultiVectorBase.hpp\
  $(srcdir)/interfaces/operator_vector/ana/extended/Thyra_ProductVectorBase.hpp\
  $(srcdir)/interfaces/operator_vector/ana/extended/Thyra_ProductVectorSpaceBase.hpp\
  $(srcdir)/interfaces/operator_vector/ana/extended/Thyra_ScaledAdjointLinearOpBaseDecl.hpp\
  $(srcdir)/interfaces/operator_vector/ana/extended/Thyra_ScaledAdjointLinearOpBase.hpp\
  $(srcdir)/interfaces/operator_vector/ana/extended/Thyra_SpmdVectorSpaceBaseDecl.hpp\
  $(srcdir)/interfaces/operator_vector/ana/extended/Thyra_SpmdVectorSpaceBase.hpp\
  $(srcdir)/interfaces/operator_vector/ana/extended/Thyra_VectorSpaceConverterBase.hpp\
  $(srcdir)/interfaces/operator_vector/ana/extended/Thyra_ZeroLinearOpBase.hpp

SUPPORT_OPERATOR_VECTOR_CLIENT_SUPPORT_H = \
  $(srcdir)/support/operator_vector/client_support/Thyra_AssertOp.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DefaultAddedLinearOpDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DefaultAddedLinearOp.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DefaultBlockedLinearOpDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DefaultBlockedLinearOp.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DefaultDiagonalLinearOpDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DefaultDiagonalLinearOp.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DefaultIdentityLinearOpDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DefaultIdentityLinearOp.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DefaultMultipliedLinearOpDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DefaultMultipliedLinearOp.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DefaultProductVectorDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DefaultProductVector.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DefaultProductVectorSpaceDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DefaultProductVectorSpace.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DefaultScaledAdjointLinearOpDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DefaultScaledAdjointLinearOp.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DefaultZeroLinearOpDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DefaultZeroLinearOp.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_describeLinearOpDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_describeLinearOp.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DetachedMultiVectorView.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_DetachedVectorView.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_LinearCombinationDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_LinearCombinationImpl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_LinearCombinationTester.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_LinearOpDefaultBaseDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_LinearOpDefaultBase.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_LinearOperatorDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_LinearOperatorImpl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_LinearOpTesterDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_LinearOpTester.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_ListedMultiVectorRandomizer.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_MultiVectorAllocator.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_MultiVectorDefaultBaseDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_MultiVectorDefaultBase.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_MultiVectorStdOpsDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_MultiVectorStdOps.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_MultiVectorStdOpsTesterDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_MultiVectorStdOpsTester.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_MultiVectorTesterDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_MultiVectorTester.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_SingleRhsLinearOpBaseDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_SingleRhsLinearOpBase.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_SingleScalarLinearOpBaseDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_SingleScalarLinearOpBase.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_SUNDIALS_Ops.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_TesterBase.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_TestingToolsDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_TestingTools.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_TestSpecifier.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_UniversalMultiVectorRandomizer.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VecOpMacros.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorDefaultBaseDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorDefaultBase.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorHandleOpsDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorHandleOpsImpl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorImpl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorOpTester.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorSpaceDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorSpaceDefaultBaseDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorSpaceDefaultBase.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorSpaceImpl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorSpaceTesterDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorSpaceTester.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorStdOpsDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorStdOps.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorStdOpsTesterDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorStdOpsTester.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorTesterDecl.hpp\
  $(srcdir)/support/operator_vector/client_support/Thyra_VectorTester.hpp

SUPPORT_OPERATOR_VECTOR_CLIENT_SUPPORT = \
  $(srcdir)/support/operator_vector/client_support/Thyra_TestingTools.cpp

SUPPORT_OPERATOR_VECTOR_ADAPTER_SUPPORT_H = \
  $(srcdir)/support/operator_vector/adapter_support/Thyra_apply_op_helper_decl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_apply_op_helper.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_DefaultClusteredSpmdProductVectorDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_DefaultClusteredSpmdProductVector.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_DefaultClusteredSpmdProductVectorSpaceDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_DefaultClusteredSpmdProductVectorSpace.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_DefaultColumnwiseMultiVectorDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_DefaultColumnwiseMultiVector.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_DefaultSpmdMultiVectorDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_DefaultSpmdMultiVector.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_DefaultSpmdVectorDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_DefaultSpmdVector.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_DefaultSpmdVectorSpaceDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_DefaultSpmdVectorSpaceFactoryDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_DefaultSpmdVectorSpaceFactory.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_DefaultSpmdVectorSpace.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_DefaultVectorMultiVectorDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_DefaultVectorMultiVector.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_EuclideanLinearOpBaseDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_EuclideanLinearOpBase.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_EuclideanScalarProdDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_EuclideanScalarProd.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_LinearOpScalarProdDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_LinearOpScalarProd.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_OperatorVectorAdapterSupportTypes.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_ScalarProdBaseDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_ScalarProdBase.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_ScalarProdVectorSpaceBaseDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_ScalarProdVectorSpaceBase.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_SingleRhsEuclideanLinearOpBaseDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_SingleRhsEuclideanLinearOpBase.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_SingleScalarEuclideanLinearOpBaseDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_SingleScalarEuclideanLinearOpBase.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_SpmdLinearOpBaseDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_SpmdLinearOpBase.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_SpmdMultiVectorBaseDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_SpmdMultiVectorBase.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_SpmdMultiVectorFileIO.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_SpmdMultiVectorSerializerDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_SpmdMultiVectorSerializer.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_SpmdVectorBaseDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_SpmdVectorBase.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_SpmdVectorSpaceDefaultBaseDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_SpmdVectorSpaceDefaultBase.hpp\
  $(srcdir)/support/operator_vector/adapter_support/Thyra_SpmdVectorSpaceUtilities.hpp

SUPPORT_OPERATOR_VECTOR_ADAPTER_SUPPORT = \
  $(srcdir)/support/operator_vector/adapter_support/Thyra_SpmdVectorSpaceUtilities.cpp

SUPPORT_OPERATOR_VECTOR_ADAPTER_SUPPORT_SERIAL_H = \
  $(srcdir)/support/operator_vector/adapter_support/serial/Thyra_DefaultSerialVectorSpaceConverterDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/serial/Thyra_DefaultSerialVectorSpaceConverter.hpp\
  $(srcdir)/support/operator_vector/adapter_support/serial/Thyra_SerialLinearOpBaseDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/serial/Thyra_SerialLinearOpBase.hpp\
  $(srcdir)/support/operator_vector/adapter_support/serial/Thyra_SerialVectorSpaceConverterBaseDecl.hpp\
  $(srcdir)/support/operator_vector/adapter_support/serial/Thyra_SerialVectorSpaceConverterBase.hpp

INTERFACES_OPERATOR_SOLVE_ANA_FUNDAMENTAL_H = \
  $(srcdir)/interfaces/operator_solve/ana/fundamental/Thyra_LinearOpSourceBase.hpp\
  $(srcdir)/interfaces/operator_solve/ana/fundamental/Thyra_LinearOpWithSolveBase.hpp\
  $(srcdir)/interfaces/operator_solve/ana/fundamental/Thyra_LinearOpWithSolveBaseDecl.hpp\
  $(srcdir)/interfaces/operator_solve/ana/fundamental/Thyra_LinearOpWithSolveFactoryBase.hpp\
  $(srcdir)/interfaces/operator_solve/ana/fundamental/Thyra_PreconditionerBase.hpp\
  $(srcdir)/interfaces/operator_solve/ana/fundamental/Thyra_PreconditionerFactoryBase.hpp\
  $(srcdir)/interfaces/operator_solve/ana/fundamental/Thyra_SolveSupportTypes.hpp

INTERFACES_OPERATOR_SOLVE_ANA_EXTENDED_H = \
  $(srcdir)/interfaces/operator_solve/ana/extended/Thyra_BlockedLinearOpWithSolveBase.hpp\
  $(srcdir)/interfaces/operator_solve/ana/extended/Thyra_InverseLinearOpBase.hpp\
  $(srcdir)/interfaces/operator_solve/ana/extended/Thyra_LinearSolverBuilderBase.hpp\
  $(srcdir)/interfaces/operator_solve/ana/extended/Thyra_PhysicallyBlockedLinearOpWithSolveBase.hpp

SUPPORT_OPERATOR_SOLVE_CLIENT_SUPPORT_H = \
  $(srcdir)/support/operator_solve/client_support/Thyra_DefaultDiagonalLinearOpWithSolve.hpp\
  $(srcdir)/support/operator_solve/client_support/Thyra_DefaultDiagonalLinearOpWithSolveDecl.hpp\
  $(srcdir)/support/operator_solve/client_support/Thyra_DefaultInverseLinearOp.hpp\
  $(srcdir)/support/operator_solve/client_support/Thyra_DefaultInverseLinearOpDecl.hpp\
  $(srcdir)/support/operator_solve/client_support/Thyra_DefaultLinearOpSource.hpp\
  $(srcdir)/support/operator_solve/client_support/Thyra_DefaultLinearOpSourceDecl.hpp\
  $(srcdir)/support/operator_solve/client_support/Thyra_DefaultPreconditioner.hpp\
  $(srcdir)/support/operator_solve/client_support/Thyra_DefaultPreconditionerDecl.hpp\
  $(srcdir)/support/operator_solve/client_support/Thyra_InverseLinearOperator.hpp\
  $(srcdir)/support/operator_solve/client_support/Thyra_LinearOpWithSolveFactoryHelpers.hpp\
  $(srcdir)/support/operator_solve/client_support/Thyra_LinearOpWithSolveTester.hpp\
  $(srcdir)/support/operator_solve/client_support/Thyra_LinearOpWithSolveTesterDecl.hpp\
  $(srcdir)/support/operator_solve/client_support/Thyra_PreconditionerFactoryHelpers.hpp\
  $(srcdir)/support/operator_solve/client_support/Thyra_SingleRhsLinearOpWithSolveBase.hpp\
  $(srcdir)/support/operator_solve/client_support/Thyra_SingleRhsLinearOpWithSolveBaseDecl.hpp\
  $(srcdir)/support/operator_solve/client_support/Thyra_SingleScalarLinearOpWithSolveBaseDecl.hpp\
  $(srcdir)/support/operator_solve/client_support/Thyra_SingleScalarLinearOpWithSolveBase.hpp

INTERFACES_NONLINEAR_MODEL_EVALUATOR_ANA_FUNDAMENTAL_H = \
  $(srcdir)/interfaces/nonlinear/model_evaluator/ana/fundamental/Thyra_ModelEvaluator.hpp\
  $(srcdir)/interfaces/nonlinear/model_evaluator/ana/fundamental/Thyra_ModelEvaluatorBase.hpp

SUPPORT_NONLINEAR_MODEL_EVALUATOR_CLIENT_SUPPORT_H = \
  $(srcdir)/support/nonlinear/model_evaluator/client_support/Thyra_DefaultEvaluationLoggerModelEvaluator.hpp\
  $(srcdir)/support/nonlinear/model_evaluator/client_support/Thyra_DefaultFinalPointCaptureModelEvaluator.hpp\
  $(srcdir)/support/nonlinear/model_evaluator/client_support/Thyra_DefaultFiniteDifferenceModelEvaluator.hpp\
  $(srcdir)/support/nonlinear/model_evaluator/client_support/Thyra_DefaultModelEvaluatorWithSolveFactory.hpp\
  $(srcdir)/support/nonlinear/model_evaluator/client_support/Thyra_DefaultMultiPeriodModelEvaluator.hpp\
  $(srcdir)/support/nonlinear/model_evaluator/client_support/Thyra_DefaultNominalBoundsOverrideModelEvaluator.hpp\
  $(srcdir)/support/nonlinear/model_evaluator/client_support/Thyra_DirectionalFiniteDiffCalculator.hpp\
  $(srcdir)/support/nonlinear/model_evaluator/client_support/Thyra_ModelEvaluatorDelegatorBase.hpp\
  $(srcdir)/support/nonlinear/model_evaluator/client_support/Thyra_ModelEvaluatorHelpers.hpp\
  $(srcdir)/support/nonlinear/model_evaluator/client_support/Thyra_SingleResidSSDAEModelEvaluator.hpp\
  $(srcdir)/support/nonlinear/model_evaluator/client_support/Thyra_StateFuncModelEvaluatorBase.hpp

INTERFACES_NONLINEAR_SOLVERS_ANA_FUNDAMENTAL_H = \
  $(srcdir)/interfaces/nonlinear/solvers/ana/fundamental/Thyra_NonlinearSolverBase.hpp

SUPPORT_NONLINEAR_SOLVERS_CLIENT_SUPPORT_H = \
  $(srcdir)/support/nonlinear/solvers/client_support/Thyra_DampenedNewtonNonlinearSolver.hpp\
  $(srcdir)/support/nonlinear/solvers/client_support/Thyra_DefaultStateEliminationModelEvaluator.hpp\
  $(srcdir)/support/nonlinear/solvers/client_support/Thyra_LinearNonlinearSolver.hpp\
  $(srcdir)/support/nonlinear/solvers/client_support/Thyra_TimeStepNewtonNonlinearSolver.hpp

#Specify locations of include files within Trilinos

include $(top_builddir)/Makefile.export.thyra

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(THYRA_INCLUDES))
else
EXPORT_INCLUDES = $(THYRA_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# --------------------------------------------------------------------
# Thyra library specifications
# --------------------------------------------------------------------

THYRA_LIB = libthyra.a

THYRA_H = \
  $(srcdir)/Thyra_ConfigDefs.hpp $(srcdir)/Thyra_Version.hpp \
  $(INTERFACES_OPERATOR_VECTOR_ANA_FUNDAMENTAL_H)\
  $(INTERFACES_OPERATOR_VECTOR_ANA_EXTENDED_H)\
  $(SUPPORT_OPERATOR_VECTOR_CLIENT_SUPPORT_H)\
  $(SUPPORT_OPERATOR_VECTOR_ADAPTER_SUPPORT_H)\
  $(SUPPORT_OPERATOR_VECTOR_ADAPTER_SUPPORT_SERIAL_H)\
  $(INTERFACES_OPERATOR_SOLVE_ANA_FUNDAMENTAL_H)\
  $(INTERFACES_OPERATOR_SOLVE_ANA_EXTENDED_H)\
  $(SUPPORT_OPERATOR_SOLVE_CLIENT_SUPPORT_H)\
  $(INTERFACES_NONLINEAR_MODEL_EVALUATOR_ANA_FUNDAMENTAL_H)\
  $(SUPPORT_NONLINEAR_MODEL_EVALUATOR_CLIENT_SUPPORT_H)\
  $(INTERFACES_NONLINEAR_SOLVERS_ANA_FUNDAMENTAL_H)\
  $(SUPPORT_NONLINEAR_SOLVERS_CLIENT_SUPPORT_H)

libthyra_a_SOURCES = $(srcdir)/Thyra_Version.cpp \
  $(SUPPORT_OPERATOR_VECTOR_CLIENT_SUPPORT) \
  $(SUPPORT_OPERATOR_VECTOR_ADAPTER_SUPPORT)

#We need the following to support separate build directories

# --------------------------------------------------------------------
# For using a special archiver
# --------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libthyra_a_AR = $(ALTERNATE_AR)

else

libthyra_a_AR = $(AR) cru

endif

# --------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# --------------------------------------------------------------------

# I have not been able to find any kind of "or" for Automake
# conditionals. JW

 libthyra_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# --------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# --------------------------------------------------------------------

lib_LIBRARIES = $(THYRA_LIB)

# --------------------------------------------------------------------
# List of all headers to install in $(includedir)
# --------------------------------------------------------------------

include_HEADERS = $(THYRA_H) 

# --------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# --------------------------------------------------------------------

# SPECIAL NOTE: ThyraConfig.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = Thyra_Config.h $(THYRA_H)

# --------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# --------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Library Target for Circular Dependencies
# ------------------------------------------------------------------------
all-libs:
	$(MAKE) libthyra.a
	cd $(top_builddir)/../rtop ; $(MAKE) all-libs
