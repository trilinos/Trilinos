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

INTREPID_LOCAL_H = \
  $(srcdir)/local/Intrepid_CellTemplates.hpp \
  $(srcdir)/local/Intrepid_RealSpaceDef.hpp \
  $(srcdir)/local/Intrepid_RealSpace.hpp \
  $(srcdir)/local/Intrepid_MultiCellDef.hpp \
  $(srcdir)/local/Intrepid_MultiCell.hpp \
  $(srcdir)/local/Intrepid_MultiCell.icpp \
  $(srcdir)/local/Intrepid_Types.hpp

INTREPID_ALL_C = \
  $(srcdir)/Intrepid_Version.cpp

# include $(top_builddir)/Makefile.export.intrepid

# The list of source files.
libintrepid_a_SOURCES = \
  $(INTREPID_ALL_C)
  
# The list of header files which are to be installed
include_HEADERS = \
  $(INTREPID_ALL_H) \
  $(INTREPID_LOCAL_H)

# All Intrepid headers are now installed in the same directory
execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = Intrepid_config.h
