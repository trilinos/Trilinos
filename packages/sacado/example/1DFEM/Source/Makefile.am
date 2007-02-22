## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#                           Sacado Package
#                 Copyright (2006) Sandia Corporation
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
# Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
# (etphipp@sandia.gov).
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# 1DFEM Sources
# ------------------------------------------------------------------------

ELEMENT_H = \
	$(srcdir)/FEApp_AbstractElement.hpp \
	$(srcdir)/FEApp_LinearElement.hpp \
	$(srcdir)/FEApp_ElementFactory.hpp

ELEMENT = \
	$(srcdir)/FEApp_LinearElement.cpp \
	$(srcdir)/FEApp_ElementFactory.cpp

QUADRATURE_H = \
	$(srcdir)/FEApp_AbstractQuadrature.hpp \
	$(srcdir)/FEApp_GaussianQuadrature2.hpp \
	$(srcdir)/FEApp_QuadratureFactory.hpp

QUADRATURE = \
	$(srcdir)/FEApp_GaussianQuadrature2.cpp \
	$(srcdir)/FEApp_QuadratureFactory.cpp

SOURCE_H = \
	$(srcdir)/FEApp_AbstractSourceFunction.hpp \
	$(srcdir)/FEApp_QuadraticSourceFunction.hpp \
	$(srcdir)/FEApp_SourceFunctionFactory.hpp \
	$(srcdir)/FEApp_SourceFunctionFactoryImpl.hpp

SOURCE = \
	$(srcdir)/FEApp_SourceFunctionFactory.cpp

PDE_H = \
	$(srcdir)/FEApp_AbstractPDE_NTBase.hpp \
	$(srcdir)/FEApp_AbstractPDE.hpp \
	$(srcdir)/FEApp_HeatNonlinearSourcePDE.hpp \
	$(srcdir)/FEApp_HeatNonlinearSourcePDEImpl.hpp \
	$(srcdir)/FEApp_BrusselatorPDE.hpp \
	$(srcdir)/FEApp_BrusselatorPDEImpl.hpp

PDE = \
	$(srcdir)/FEApp_HeatNonlinearSourcePDE.cpp \
	$(srcdir)/FEApp_BrusselatorPDE.cpp

BC_H = \
	$(srcdir)/FEApp_AbstractBC.hpp \
	$(srcdir)/FEApp_ConstantDirichletBC.hpp

BC = \
	$(srcdir)/FEApp_ConstantDirichletBC.cpp

PROBLEM_H = \
	$(srcdir)/FEApp_TemplateTypes.hpp \
	$(srcdir)/FEApp_AbstractProblem.hpp \
	$(srcdir)/FEApp_BrusselatorProblem.hpp \
	$(srcdir)/FEApp_HeatNonlinearSourceProblem.hpp \
	$(srcdir)/FEApp_ProblemFactory.hpp

PROBLEM = \
	$(srcdir)/FEApp_BrusselatorProblem.cpp \
	$(srcdir)/FEApp_HeatNonlinearSourceProblem.cpp \
	$(srcdir)/FEApp_ProblemFactory.cpp

MESH_H = \
	$(srcdir)/FEApp_Mesh.hpp

MESH = \
	$(srcdir)/FEApp_Mesh.cpp

FILL_H = \
	$(srcdir)/FEApp_AbstractInitPostOp.hpp \
	$(srcdir)/FEApp_InitPostOps.hpp \
	$(srcdir)/FEApp_GlobalFill.hpp \
	$(srcdir)/FEApp_GlobalFillImpl.hpp

FILL = \
	$(srcdir)/FEApp_InitPostOps.cpp \
	$(srcdir)/FEApp_GlobalFill.cpp

DISC_H = \
	$(srcdir)/FEApp_AbstractDiscretization.hpp \
	$(srcdir)/FEApp_CZeroDiscretization.hpp \
	$(srcdir)/FEApp_DiscretizationFactory.hpp

DISC = \
	$(srcdir)/FEApp_CZeroDiscretization.cpp \
	$(srcdir)/FEApp_DiscretizationFactory.cpp

APP_H = \
	$(srcdir)/FEApp_Application.hpp \
	$(srcdir)/FEApp_ModelEvaluator.hpp

APP = \
	$(srcdir)/FEApp_Application.cpp \
	$(srcdir)/FEApp_ModelEvaluator.cpp

ONEDFEM_H = \
	$(ELEMENT_H) \
	$(QUADRATURE_H) \
	$(SOURCE_H) \
	$(PDE_H) \
	$(BC_H) \
	$(PROBLEM_H) \
	$(MESH_H) \
	$(FILL_H) \
	$(DISC_H) \
	$(APP_H)

ONEDFEM = \
	$(ELEMENT) \
	$(QUADRATURE) \
	$(SOURCE) \
	$(PDE) \
	$(BC) \
	$(PROBLEM) \
	$(MESH) \
	$(FILL) \
	$(DISC) \
	$(APP)

# ------------------------------------------------------------------------
# 1DFEM library specifications
# ------------------------------------------------------------------------

lib_LIBRARIES = lib1DFEM.a

lib1DFEM_a_SOURCES = $(ONEDFEM) $(ONEDFEM_H)

include_HEADERS = $(ONEDFEM_H) 

#We need the following to support separate build directories
AM_CPPFLAGS=\
	-I$(top_srcdir)/src \
	-I$(top_srcdir)/src/mpl \
	-I$(top_builddir)/src \
	$(LOCA_INCLUDES)

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

lib1DFEM_a_AR = $(ALTERNATE_AR)

else

lib1DFEM_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

lib1DFEM_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# Include NOX export Makefile for correct paths
include $(top_builddir)/../nox/Makefile.export.nox

