## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#                           Sacado Package
#                 Copyright (2004) Sandia Corporation
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
# Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
# (etphipp@sandia.gov).
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# 1DFEM Sources
# ------------------------------------------------------------------------

ONEDFEM_H = \
	$(srcdir)/FEApp_AbstractElement.hpp \
	$(srcdir)/FEApp_LinearElement.hpp \
	$(srcdir)/FEApp_AbstractQuadrature.hpp \
	$(srcdir)/FEApp_GaussianQuadrature2.hpp \
	$(srcdir)/FEApp_AbstractSourceFunction.hpp \
	$(srcdir)/FEApp_QuadraticSourceFunction.hpp \
	$(srcdir)/FEApp_AbstractPDE.hpp \
	$(srcdir)/FEApp_HeatNonlinearSourcePDE.hpp \
	$(srcdir)/FEApp_HeatNonlinearSourcePDEImpl.hpp \
	$(srcdir)/FEApp_AbstractBC.hpp \
	$(srcdir)/FEApp_ConstantDirichletBC.hpp \
	$(srcdir)/FEApp_Mesh.hpp \
	$(srcdir)/FEApp_AbstractInitPostOp.hpp \
	$(srcdir)/FEApp_InitPostOps.hpp \
	$(srcdir)/FEApp_GlobalFill.hpp \
	$(srcdir)/FEApp_GlobalFillImpl.hpp \
	$(srcdir)/FEApp_CZeroDiscretization.hpp \
	$(srcdir)/FEApp_Application.hpp

ONEDFEM = \
	$(srcdir)/FEApp_LinearElement.cpp \
	$(srcdir)/FEApp_GaussianQuadrature2.cpp \
	$(srcdir)/FEApp_ConstantDirichletBC.cpp \
	$(srcdir)/FEApp_Mesh.cpp \
	$(srcdir)/FEApp_InitPostOps.cpp \
	$(srcdir)/FEApp_CZeroDiscretization.cpp \
	$(srcdir)/FEApp_Application.cpp

# ------------------------------------------------------------------------
# 1DFEM library specifications
# ------------------------------------------------------------------------

lib_LIBRARIES = lib1DFEM.a

lib1DFEM_a_SOURCES = $(ONEDFEM) $(ONEDFEM_H)

include_HEADERS = $(ONEDFEM_H) 

#We need the following to support separate build directories
AM_CPPFLAGS=\
	-I$(top_srcdir)/src \
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

