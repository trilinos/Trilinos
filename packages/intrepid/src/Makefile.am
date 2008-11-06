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

INTREPID_COMMON_H = \
  $(srcdir)/common/Intrepid_FieldContainerDef.hpp \
  $(srcdir)/common/Intrepid_FieldContainer.hpp \
  $(srcdir)/common/Intrepid_RealSpaceDef.hpp \
  $(srcdir)/common/Intrepid_RealSpace.hpp \
  $(srcdir)/common/Intrepid_RealSpaceToolsDef.hpp \
  $(srcdir)/common/Intrepid_RealSpaceTools.hpp \
  $(srcdir)/common/Intrepid_Types.hpp \
  $(srcdir)/common/Intrepid_Utils.hpp

INTREPID_LOCAL_H =  \
  $(srcdir)/local/basis/Intrepid_Basis.hpp \
  $(srcdir)/local/basis/Intrepid_DefaultBasisFactoryDef.hpp \
  $(srcdir)/local/basis/Intrepid_DefaultBasisFactory.hpp \
  $(srcdir)/local/basis/Intrepid_Tabulate.hpp \
  $(srcdir)/local/basis/Intrepid_F0_QUAD_I1_FEM_DEFAULTDef.hpp \
  $(srcdir)/local/basis/Intrepid_F0_QUAD_I1_FEM_DEFAULT.hpp \
  $(srcdir)/local/basis/Intrepid_F0_TRI_C1_FEM_DEFAULTDef.hpp \
  $(srcdir)/local/basis/Intrepid_F0_TRI_C1_FEM_DEFAULT.hpp \
  $(srcdir)/local/basis/Intrepid_F0_TRI_C2_FEM_DEFAULTDef.hpp \
  $(srcdir)/local/basis/Intrepid_F0_TRI_C2_FEM_DEFAULT.hpp \
  $(srcdir)/local/Intrepid_CellDef.hpp \
  $(srcdir)/local/Intrepid_Cell.hpp \
  $(srcdir)/local/Intrepid_CellTemplates.hpp \
  $(srcdir)/local/Intrepid_CubatureDirectDef.hpp \
  $(srcdir)/local/Intrepid_CubatureDirect.hpp \
  $(srcdir)/local/Intrepid_Cubature.hpp \
  $(srcdir)/local/Intrepid_CubatureGenSparseDef.hpp \
  $(srcdir)/local/Intrepid_CubatureGenSparse.hpp \
  $(srcdir)/local/Intrepid_CubatureSparseDef.hpp \
  $(srcdir)/local/Intrepid_CubatureSparse.hpp \
  $(srcdir)/local/Intrepid_CubatureSparseHelper.hpp \
  $(srcdir)/local/Intrepid_CubatureTensorDef.hpp \
  $(srcdir)/local/Intrepid_CubatureTensor.hpp \
  $(srcdir)/local/Intrepid_DefaultCubatureFactory.hpp \
  $(srcdir)/local/Intrepid_DefaultCubatureFactoryDef.hpp \
  $(srcdir)/local/Intrepid_DefaultFieldFactory.hpp \
  $(srcdir)/local/Intrepid_DefaultFieldFactoryDef.hpp \
  $(srcdir)/local/Intrepid_LocalField.hpp \
  $(srcdir)/local/Intrepid_LocalForm0Def.hpp \
  $(srcdir)/local/Intrepid_LocalForm0.hpp \
  $(srcdir)/local/Intrepid_MultiCellDef.hpp \
  $(srcdir)/local/Intrepid_MultiCell.hpp \
  $(srcdir)/local/Intrepid_MultiCell.icpp \
  $(srcdir)/local/basis/Intrepid_F0_TRI_C1_FEM_FIAT.hpp \
  $(srcdir)/local/basis/Intrepid_F0_TRI_C1_FEM_FIATDef.hpp \
  $(srcdir)/local/basis/Intrepid_F0_TRI_C2_FEM_FIAT.hpp \
  $(srcdir)/local/basis/Intrepid_F0_TRI_C2_FEM_FIATDef.hpp \
  $(srcdir)/local/basis/Intrepid_F0_TRI_C3_FEM_FIAT.hpp \
  $(srcdir)/local/basis/Intrepid_F0_TRI_C3_FEM_FIATDef.hpp \
  $(srcdir)/local/basis/Intrepid_F0_TET_C1_FEM_FIAT.hpp \
  $(srcdir)/local/basis/Intrepid_F0_TET_C1_FEM_FIATDef.hpp \
  $(srcdir)/local/basis/Intrepid_F0_TET_C2_FEM_FIAT.hpp \
  $(srcdir)/local/basis/Intrepid_F0_TET_C2_FEM_FIATDef.hpp \
  $(srcdir)/local/basis/Intrepid_F0_TET_C3_FEM_FIAT.hpp \
  $(srcdir)/local/basis/Intrepid_F0_TET_C3_FEM_FIATDef.hpp \
  $(srcdir)/local/basis/Intrepid_F2_TRI_I1_FEM_FIAT.hpp \
  $(srcdir)/local/basis/Intrepid_F2_TRI_I1_FEM_FIATDef.hpp \
  $(srcdir)/local/basis/Intrepid_F2_TRI_I2_FEM_FIAT.hpp \
  $(srcdir)/local/basis/Intrepid_F2_TRI_I2_FEM_FIATDef.hpp \
  $(srcdir)/local/basis/Intrepid_F2_TRI_I3_FEM_FIAT.hpp \
  $(srcdir)/local/basis/Intrepid_F2_TRI_I3_FEM_FIATDef.hpp \
  $(srcdir)/local/basis/Intrepid_F2_TET_I1_FEM_FIAT.hpp \
  $(srcdir)/local/basis/Intrepid_F2_TET_I1_FEM_FIATDef.hpp \
  $(srcdir)/local/basis/Intrepid_F2_TET_I2_FEM_FIAT.hpp \
  $(srcdir)/local/basis/Intrepid_F2_TET_I2_FEM_FIATDef.hpp \
  $(srcdir)/local/basis/Intrepid_F2_TET_I3_FEM_FIAT.hpp \
  $(srcdir)/local/basis/Intrepid_F2_TET_I3_FEM_FIATDef.hpp \
  $(srcdir)/local/basis/Intrepid_F1_TRI_I1_FEM_FIAT.hpp \
  $(srcdir)/local/basis/Intrepid_F1_TRI_I1_FEM_FIATDef.hpp \
  $(srcdir)/local/basis/Intrepid_F1_TRI_I2_FEM_FIAT.hpp \
  $(srcdir)/local/basis/Intrepid_F1_TRI_I2_FEM_FIATDef.hpp \
  $(srcdir)/local/basis/Intrepid_F1_TRI_I3_FEM_FIAT.hpp \
  $(srcdir)/local/basis/Intrepid_F1_TRI_I3_FEM_FIATDef.hpp \
  $(srcdir)/local/basis/Intrepid_F1_TET_I1_FEM_FIAT.hpp \
  $(srcdir)/local/basis/Intrepid_F1_TET_I1_FEM_FIATDef.hpp \
  $(srcdir)/local/basis/Intrepid_F1_TET_I2_FEM_FIAT.hpp \
  $(srcdir)/local/basis/Intrepid_F1_TET_I2_FEM_FIATDef.hpp \
  $(srcdir)/local/basis/Intrepid_F1_TET_I3_FEM_FIAT.hpp \
  $(srcdir)/local/basis/Intrepid_F1_TET_I3_FEM_FIATDef.hpp \
  $(srcdir)/local/basis/FIAT/FIAT_AffineMap.hpp \
  $(srcdir)/local/basis/FIAT/FIAT_AffineMapDef.hpp \
  $(srcdir)/local/basis/FIAT/FIAT_ReferenceElement.hpp \
  $(srcdir)/local/basis/FIAT/FIAT_ReferenceElementDef.hpp \
  $(srcdir)/local/basis/FIAT/FIAT_Line.hpp \
  $(srcdir)/local/basis/FIAT/FIAT_LineDef.hpp \
  $(srcdir)/local/basis/FIAT/FIAT_Line.hpp \
  $(srcdir)/local/basis/FIAT/FIAT_LineDef.hpp \
  $(srcdir)/local/basis/FIAT/FIAT_DefaultLine.hpp \
  $(srcdir)/local/basis/FIAT/FIAT_DefaultLineDef.hpp \
  $(srcdir)/local/basis/FIAT/Lagrange1d.hpp \
  $(srcdir)/local/basis/FIAT/Lagrange1dDef.hpp \
  $(srcdir)/local/basis/Intrepid_F0_QUAD_DD.hpp \
  $(srcdir)/local/basis/Intrepid_F0_QUAD_DDDef.hpp \
  $(srcdir)/local/basis/Intrepid_F0_HEX_DD.hpp \
  $(srcdir)/local/basis/Intrepid_F0_HEX_DDDef.hpp


INTREPID_ALL_C = \
  $(srcdir)/Intrepid_Version.cpp

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
  $(INTREPID_COMMON_H) \
  $(INTREPID_LOCAL_H)

# All Intrepid headers are now installed in the same directory
execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = Intrepid_config.h









