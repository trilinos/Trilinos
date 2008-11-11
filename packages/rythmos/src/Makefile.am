#@HEADER
# ***********************************************************************
#
#                           Rythmos Package
#                 Copyright (2006) Sandia Corporation
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
# Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
#
# ***********************************************************************
#@HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package rythmos subdirectory src built successfully."
	@echo ""

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

CORE_H = \
  $(srcdir)/Rythmos_AdjointModelEvaluator.hpp\
  $(srcdir)/Rythmos_BackwardEulerStepper.hpp\
  $(srcdir)/Rythmos_BreakPointInformerBase.hpp\
  $(srcdir)/Rythmos_CompositeIntegratorObserver.hpp\
  $(srcdir)/Rythmos_ConfigDefs.h\
  $(srcdir)/Rythmos_DataStore.hpp\
  $(srcdir)/Rythmos_DefaultIntegrator.hpp\
	$(srcdir)/Rythmos_DiagonalImplicitRKModelEvaluator.hpp\
	$(srcdir)/Rythmos_ErrorEstimateBase.hpp\
	$(srcdir)/Rythmos_ErrorEstimatorBase.hpp\
  $(srcdir)/Rythmos_ErrWtVecCalcAcceptingStepControlStrategyBase.hpp\
  $(srcdir)/Rythmos_ErrWtVecCalcBase.hpp\
  $(srcdir)/Rythmos_ExplicitRKStepper.hpp\
  $(srcdir)/Rythmos_ExplicitTaylorPolynomialStepper.hpp\
  $(srcdir)/Rythmos_ForwardEulerStepper.hpp\
  $(srcdir)/Rythmos_ForwardResponseSensitivityComputer.hpp\
  $(srcdir)/Rythmos_ForwardResponseSensitivityComputerObserver.hpp\
  $(srcdir)/Rythmos_ForwardSensitivityIntegratorAsModelEvaluator.hpp\
  $(srcdir)/Rythmos_ForwardSensitivityModelEvaluator.hpp\
  $(srcdir)/Rythmos_ForwardSensitivityStepper.hpp\
  $(srcdir)/Rythmos_HermiteInterpolator.hpp\
  $(srcdir)/Rythmos_ImplicitBDFStepper.hpp\
  $(srcdir)/Rythmos_ImplicitBDFStepperDecl.hpp\
  $(srcdir)/Rythmos_ImplicitBDFStepperErrWtVecCalc.hpp\
  $(srcdir)/Rythmos_ImplicitBDFStepperStepControl.hpp\
  $(srcdir)/Rythmos_ImplicitBDFStepperStepControlDecl.hpp\
  $(srcdir)/Rythmos_ImplicitRKModelEvaluator.hpp\
  $(srcdir)/Rythmos_ImplicitRKStepper.hpp\
  $(srcdir)/Rythmos_IntegrationControlStrategyBase.hpp\
  $(srcdir)/Rythmos_IntegrationObserverBase.hpp\
  $(srcdir)/Rythmos_IntegratorBase.hpp\
  $(srcdir)/Rythmos_InterpolationBuffer.hpp\
  $(srcdir)/Rythmos_InterpolationBufferAppenderBase.hpp\
  $(srcdir)/Rythmos_InterpolationBufferBase.hpp\
  $(srcdir)/Rythmos_InterpolationBufferHelpers.hpp\
  $(srcdir)/Rythmos_InterpolatorBase.hpp\
  $(srcdir)/Rythmos_LinearInterpolator.hpp\
  $(srcdir)/Rythmos_PointwiseInterpolationBufferAppender.hpp\
  $(srcdir)/Rythmos_QuadratureBase.hpp\
  $(srcdir)/Rythmos_RKButcherTableau.hpp\
  $(srcdir)/Rythmos_ResponseAndFwdSensPoint.hpp\
  $(srcdir)/Rythmos_SimpleIntegrationControlStrategy.hpp\
  $(srcdir)/Rythmos_SingleResidualModelEvaluator.hpp\
  $(srcdir)/Rythmos_SingleResidualModelEvaluatorBase.hpp\
  $(srcdir)/Rythmos_SmartInterpolationBufferAppender.hpp\
  $(srcdir)/Rythmos_SolverAcceptingStepperBase.hpp\
  $(srcdir)/Rythmos_StateAndForwardSensitivityModelEvaluator.hpp\
  $(srcdir)/Rythmos_StepControlInfo.hpp\
  $(srcdir)/Rythmos_StepControlStrategyAcceptingStepperBase.hpp\
  $(srcdir)/Rythmos_StepControlStrategyBase.hpp\
  $(srcdir)/Rythmos_StepperAsModelEvaluator.hpp\
  $(srcdir)/Rythmos_StepperBase.hpp\
  $(srcdir)/Rythmos_StepperHelpers.hpp\
  $(srcdir)/Rythmos_StepperSupportTypes.hpp\
  $(srcdir)/Rythmos_TimeDiscretizedBackwardEulerModelEvaluator.hpp\
  $(srcdir)/Rythmos_TimeRange.hpp\
  $(srcdir)/Rythmos_TimeStepNonlinearSolver.hpp\
  $(srcdir)/Rythmos_Types.hpp\
  $(srcdir)/Rythmos_extractStateAndSens.hpp\
  $(srcdir)/Rythmos_Version.h


CORE = \
  $(srcdir)/Rythmos_Version.cpp\
	$(srcdir)/Rythmos_RKButcherTableau.cpp

UTIL_H =

UTIL = 

ENUM_H = 

UNIPROC_H = 

UNIPROC = 

SERIAL_H = 

SERIAL = 

EXTRA_MPI_H =

EXTRA_MPI = 

if HAVE_MPI 
MPI_H = \
  $(EXTRA_MPI_H)
MPI = \
  $(EXTRA_MPI)

else

MPI_H =
MPI =

endif

USERBASE_H = 


LANG_WRAPS_H = 

LANG_WRAPS = 

DIST_H = 

DIST = 

COMP_H =

COMP = 

#Specify locations of include files within Trilinos

include $(top_builddir)/Makefile.export.rythmos

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(RYTHMOS_INCLUDES))
else
EXPORT_INCLUDES = $(RYTHMOS_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# ------------------------------------------------------------------------
# Rythmos library specifications
# ------------------------------------------------------------------------

RYTHMOS_LIB = librythmos.a

RYTHMOS_H = \
  $(CORE_H) \
  $(UTIL_H) \
  $(ENUM_H) \
  $(UNIPROC_H) \
  $(SERIAL_H) \
  $(USERBASE_H) \
  $(DIST_H) \
  $(COMP_H) \
  $(MPI_H)

librythmos_a_SOURCES = \
  $(CORE) \
  $(UTIL) \
  $(UNIPROC) \
  $(SERIAL) \
  $(DIST) \
  $(COMP) \
  $(MPI)

EXTRA_librythmos_a_SOURCES = \
  $(EXTRA_MPI) \
  $(EXTRA_MPI_H)

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

librythmos_a_AR = $(ALTERNATE_AR)
else

librythmos_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

librythmos_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(RYTHMOS_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(RYTHMOS_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# All Trilinos headers are now installed in the same directory
execincludedir = $(includedir)

nodist_execinclude_HEADERS = Rythmos_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 
