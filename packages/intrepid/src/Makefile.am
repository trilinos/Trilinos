# @HEADER
# ************************************************************************
#
#                           Intrepid Package
#                 Copyright (2007) Sandia Corporation
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
# Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
#                    Denis Ridzal (dridzal@sandia.gov).
#
# ************************************************************************
# @HEADER

all-local:
	@echo ""
	@echo "Trilinos package Intrepid subdirectory src built successfully."
	@echo ""


# Subdirectories to build
SUBDIRS =

# Things to build
lib_LIBRARIES = libintrepid.a

INTREPID_ALL_H = \
  $(srcdir)/Intrepid_ConfigDefs.hpp

INTREPID_SHARED_H = \
  $(srcdir)/Shared/Intrepid_ArrayToolsDef.hpp \
  $(srcdir)/Shared/Intrepid_ArrayTools.hpp \
  $(srcdir)/Shared/Intrepid_FieldContainerDef.hpp \
  $(srcdir)/Shared/Intrepid_FieldContainer.hpp \
  $(srcdir)/Shared/Intrepid_PolylibDef.hpp \
  $(srcdir)/Shared/Intrepid_Polylib.hpp \
  $(srcdir)/Shared/Intrepid_RealSpaceToolsDef.hpp \
  $(srcdir)/Shared/Intrepid_RealSpaceTools.hpp \
  $(srcdir)/Shared/Intrepid_Types.hpp \
  $(srcdir)/Shared/Intrepid_Utils.hpp

INTREPID_DISCRETIZATION_H = \
  $(srcdir)/Discretization/Integration/Intrepid_Cubature.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_CubatureDirectDef.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_CubatureDirect.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_CubatureDirectLineGaussDef.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_CubatureDirectLineGauss.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_CubatureDirectTetDefaultDef.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_CubatureDirectTetDefault.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_CubatureDirectTriDefaultDef.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_CubatureDirectTriDefault.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_CubatureGenSparseDef.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_CubatureGenSparse.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_CubatureSparseDef.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_CubatureSparse.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_CubatureSparseHelper.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_CubatureTensorDef.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_CubatureTensor.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_CubaturePolylib.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_CubaturePolylibDef.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_DefaultCubatureFactoryDef.hpp \
  $(srcdir)/Discretization/Integration/Intrepid_DefaultCubatureFactory.hpp \
  $(srcdir)/Discretization/FunctionSpaceTools/Intrepid_FunctionSpaceToolsDef.hpp \
  $(srcdir)/Discretization/FunctionSpaceTools/Intrepid_FunctionSpaceTools.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_BasisDef.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_Basis.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HGRAD_HEX_C1_FEM.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HGRAD_HEX_C1_FEMDef.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HGRAD_QUAD_C1_FEM.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HGRAD_QUAD_C1_FEMDef.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HGRAD_TET_C1_FEM.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HGRAD_TET_C1_FEMDef.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HGRAD_TRI_C1_FEM.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HGRAD_TRI_C1_FEMDef.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HGRAD_WEDGE_C1_FEM.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HGRAD_WEDGE_C1_FEMDef.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HCURL_HEX_I1_FEM.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HCURL_HEX_I1_FEMDef.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HCURL_TET_I1_FEM.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HCURL_TET_I1_FEMDef.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HCURL_WEDGE_I1_FEM.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HCURL_WEDGE_I1_FEMDef.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HDIV_HEX_I1_FEM.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HDIV_HEX_I1_FEMDef.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HDIV_TET_I1_FEM.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HDIV_TET_I1_FEMDef.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HDIV_WEDGE_I1_FEM.hpp \
  $(srcdir)/Discretization/Basis/Intrepid_HDIV_WEDGE_I1_FEMDef.hpp 

INTREPID_CELL_H = \
  $(srcdir)/Cell/Intrepid_CellTools.hpp \
  $(srcdir)/Cell/Intrepid_CellToolsDef.hpp


INTREPID_ALL_C = \
  $(srcdir)/Intrepid_Version.cpp \
  $(srcdir)/Shared/Intrepid_Utils.cpp

# ------------------------------------------------------------------------
# Includes for required and optional library dependencies
# ------------------------------------------------------------------------
# List locations of include files within Trilinos but
# outside of Intrepid.
include $(top_builddir)/Makefile.export.intrepid

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(INTREPID_INCLUDES))
else
EXPORT_INCLUDES = $(INTREPID_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# The list of source files.
libintrepid_a_SOURCES = \
  $(INTREPID_ALL_C)

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------
# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.
if USE_ALTERNATE_AR
libintrepid_a_AR = $(ALTERNATE_AR)
else
libintrepid_a_AR = $(AR) cru
endif
  
# The list of header files which are to be installed
include_HEADERS = \
  $(INTREPID_ALL_H) \
  $(INTREPID_SHARED_H) \
  $(INTREPID_DISCRETIZATION_H) \
  $(INTREPID_CELL_H)

# All Intrepid headers are now installed in the same directory
execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = Intrepid_config.h









