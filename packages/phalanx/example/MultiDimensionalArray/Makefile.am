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

AM_CPPFLAGS = -I$(srcdir)/evaluators $(EXAMPLE_INCLUDES)

EXEEXT = .exe

if BUILD_EXAMPLES
noinst_PROGRAMS = example_multidimensional_array
endif

example_multidimensional_array_SOURCES = \
	$(srcdir)/Example_MultiDimensionalArray.cpp \
	$(srcdir)/Traits.hpp \
	$(srcdir)/Traits.cpp \
	$(srcdir)/FactoryTraits.hpp \
	$(srcdir)/Cell.hpp \
	$(srcdir)/Cell.cpp \
	$(srcdir)/Dimension.hpp \
	$(srcdir)/Dimension.cpp \
	$(srcdir)/Workset.hpp \
	$(srcdir)/evaluators/Evaluator_Constant.hpp \
	$(srcdir)/evaluators/Evaluator_Constant_Def.hpp \
	$(srcdir)/evaluators/Evaluator_Density.hpp \
	$(srcdir)/evaluators/Evaluator_Density_Def.hpp \
	$(srcdir)/evaluators/Evaluator_EnergyFlux_Fourier.hpp \
	$(srcdir)/evaluators/Evaluator_EnergyFlux_Fourier_Def.hpp \
	$(srcdir)/evaluators/Evaluator_FEInterpolation.hpp \
	$(srcdir)/evaluators/Evaluator_FEInterpolation_Def.hpp \
	$(srcdir)/evaluators/Evaluator_NonlinearSource.hpp \
	$(srcdir)/evaluators/Evaluator_NonlinearSource_Def.hpp

example_multidimensional_array_DEPENDENCIES = \
	$(top_builddir)/src/libphalanx.a

example_multidimensional_array_LDADD = $(EXAMPLE_LIBS)
