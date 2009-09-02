#@HEADER
# ************************************************************************
# 
#        Phalanx: A Partial Differential Equation Field Evaluation 
#       Kernel for Flexible Management of Complex Dependency Chains
#                  Copyright (2008) Sandia Corporation
# 
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
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
# Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
# National Laboratories.
# 
# ************************************************************************
#@HEADER

include $(top_builddir)/Makefile.export.phalanx

# The following line helps the test harness recover from build errors.

all-local:

if USING_GNUMAKE
EXAMPLE_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(PHALANX_LIBS))
EXAMPLE_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(PHALANX_INCLUDES))
else
EXAMPLE_LIBS     = $(PHALANX_LIBS)
EXAMPLE_INCLUDES = $(PHALANX_INCLUDES)
endif

AM_CPPFLAGS = -I$(srcdir)/evaluators -I$(srcdir)/../FEM_Nonlinear $(EXAMPLE_INCLUDES)

EXEEXT = .exe

if BUILD_PHALANX_BELOS
if BUILD_PHALANX_IFPACK
if BUILD_PHALANX_EPETRA
if BUILD_EXAMPLES
noinst_PROGRAMS = example_fem_nonlinear_vbr example_epetra_vbr_test
endif
endif
endif
endif

example_fem_nonlinear_vbr_SOURCES = \
	$(srcdir)/Example_FEM_Nonlinear_VBR.cpp \
	$(srcdir)/LinearObjectFactoryVBR.hpp \
	$(srcdir)/LinearObjectFactoryVBR.cpp \
	$(srcdir)/../FEM_Nonlinear/Traits.hpp \
	$(srcdir)/../FEM_Nonlinear/Traits.cpp \
	$(srcdir)/../FEM_Nonlinear/FactoryTraits.hpp \
	$(srcdir)/../FEM_Nonlinear/Element_Linear2D.hpp \
	$(srcdir)/../FEM_Nonlinear/Element_Linear2D.cpp \
	$(srcdir)/../FEM_Nonlinear/ExplicitTemplateInstantiation.hpp \
	$(srcdir)/../FEM_Nonlinear/Dimension.hpp \
	$(srcdir)/../FEM_Nonlinear/Dimension.cpp \
	$(srcdir)/../FEM_Nonlinear/MeshBuilder.hpp \
	$(srcdir)/../FEM_Nonlinear/MeshBuilder.cpp \
	$(srcdir)/../FEM_Nonlinear/Workset.hpp \
	$(srcdir)/evaluators/Evaluator_Constant.hpp \
	$(srcdir)/evaluators/Evaluator_Constant_Def.hpp \
	$(srcdir)/evaluators/Evaluator_Constant.cpp \
	$(srcdir)/evaluators/Evaluator_Equations.hpp \
	$(srcdir)/evaluators/Evaluator_Equations_Def.hpp \
	$(srcdir)/evaluators/Evaluator_Equations.cpp \
	$(srcdir)/evaluators/Evaluator_FEInterpolation.hpp \
	$(srcdir)/evaluators/Evaluator_FEInterpolation_Def.hpp \
	$(srcdir)/evaluators/Evaluator_FEInterpolation.cpp \
	$(srcdir)/evaluators/Evaluator_GatherSolution.hpp \
	$(srcdir)/evaluators/Evaluator_GatherSolution_Def.hpp \
	$(srcdir)/evaluators/Evaluator_GatherSolution.cpp \
	$(srcdir)/evaluators/Evaluator_ScatterResidual.hpp \
	$(srcdir)/evaluators/Evaluator_ScatterResidual_Def.hpp \
	$(srcdir)/evaluators/Evaluator_ScatterResidual.cpp

example_fem_nonlinear_vbr_DEPENDENCIES = \
	$(top_builddir)/src/libphalanx.a

example_fem_nonlinear_vbr_LDADD = $(EXAMPLE_LIBS)

example_epetra_vbr_test_SOURCES = \
	$(srcdir)/auxiliary_tests/Example_Epetra_VBR_Test.cpp \
	$(srcdir)/../FEM_Nonlinear/Traits.hpp \
	$(srcdir)/../FEM_Nonlinear/Traits.cpp \
	$(srcdir)/../FEM_Nonlinear/FactoryTraits.hpp \
	$(srcdir)/../FEM_Nonlinear/Element_Linear2D.hpp \
	$(srcdir)/../FEM_Nonlinear/Element_Linear2D.cpp \
	$(srcdir)/../FEM_Nonlinear/ExplicitTemplateInstantiation.hpp \
	$(srcdir)/../FEM_Nonlinear/Dimension.hpp \
	$(srcdir)/../FEM_Nonlinear/Dimension.cpp \
	$(srcdir)/../FEM_Nonlinear/Workset.hpp \
	$(srcdir)/evaluators/Evaluator_Constant.hpp \
	$(srcdir)/evaluators/Evaluator_Constant_Def.hpp \
	$(srcdir)/evaluators/Evaluator_Constant.cpp \
	$(srcdir)/evaluators/Evaluator_Equations.hpp \
	$(srcdir)/evaluators/Evaluator_Equations_Def.hpp \
	$(srcdir)/evaluators/Evaluator_Equations.cpp \
	$(srcdir)/evaluators/Evaluator_FEInterpolation.hpp \
	$(srcdir)/evaluators/Evaluator_FEInterpolation_Def.hpp \
	$(srcdir)/evaluators/Evaluator_FEInterpolation.cpp \
	$(srcdir)/evaluators/Evaluator_GatherSolution.hpp \
	$(srcdir)/evaluators/Evaluator_GatherSolution_Def.hpp \
	$(srcdir)/evaluators/Evaluator_GatherSolution.cpp \
	$(srcdir)/evaluators/Evaluator_ScatterResidual.hpp \
	$(srcdir)/evaluators/Evaluator_ScatterResidual_Def.hpp \
	$(srcdir)/evaluators/Evaluator_ScatterResidual.cpp

example_epetra_vbr_test_DEPENDENCIES = \
	$(top_builddir)/src/libphalanx.a

example_epetra_vbr_test_LDADD = $(EXAMPLE_LIBS)
