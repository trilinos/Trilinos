# @HEADER
# ************************************************************************
#
#           Galeri: Finite Element and Matrix Generation Package
#                 Copyright (2006) ETHZ/Sandia Corporation
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
# Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
#
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package galeri subdirectory src built successfully."
	@echo ""

CORE_H = $(srcdir)/Galeri_Version.h                    \
	 $(srcdir)/Galeri_ConfigDefs.h                 \
	 $(srcdir)/Galeri_Maps.h                       \
	 $(srcdir)/Galeri_CrsMatrices.h                \
	 $(srcdir)/Galeri_VbrMatrices.h                \
         $(srcdir)/Galeri_Exception.h                  \
	 $(srcdir)/Galeri_FiniteElements.h             \
	 $(srcdir)/Galeri_Utils.h                      \
	 $(srcdir)/Galeri_ReadHB.h                     \
	 $(srcdir)/iohb.h                              \
	 $(srcdir)/Maps/Galeri_Linear.h                \
	 $(srcdir)/Maps/Galeri_Interlaced.h            \
	 $(srcdir)/Maps/Galeri_Random.h                \
         $(srcdir)/Maps/Galeri_Cartesian2D.h           \
         $(srcdir)/Maps/Galeri_Cartesian3D.h           \
	 $(srcdir)/CrsMatrices/Galeri_BigCross2D.h        \
	 $(srcdir)/CrsMatrices/Galeri_BigStar2D.h         \
	 $(srcdir)/CrsMatrices/Galeri_Cauchy.h            \
	 $(srcdir)/CrsMatrices/Galeri_KMS.h               \
	 $(srcdir)/CrsMatrices/Galeri_Ris.h               \
	 $(srcdir)/CrsMatrices/Galeri_Cross2D.h           \
	 $(srcdir)/CrsMatrices/Galeri_Laplace1DNeumann.h  \
	 $(srcdir)/CrsMatrices/Galeri_Star2D.h            \
	 $(srcdir)/CrsMatrices/Galeri_Cross3D.h           \
	 $(srcdir)/CrsMatrices/Galeri_Lehmer.h            \
	 $(srcdir)/CrsMatrices/Galeri_Stretched2D.h       \
	 $(srcdir)/CrsMatrices/Galeri_Diag.h              \
	 $(srcdir)/CrsMatrices/Galeri_Minij.h             \
	 $(srcdir)/CrsMatrices/Galeri_Tridiag.h           \
	 $(srcdir)/CrsMatrices/Galeri_Fielder.h           \
	 $(srcdir)/CrsMatrices/Galeri_Ones.h              \
	 $(srcdir)/CrsMatrices/Galeri_UniFlow2D.h         \
	 $(srcdir)/CrsMatrices/Galeri_Hanowa.h            \
	 $(srcdir)/CrsMatrices/Galeri_Parter.h            \
	 $(srcdir)/CrsMatrices/Galeri_Vander.h            \
	 $(srcdir)/CrsMatrices/Galeri_Hilbert.h           \
	 $(srcdir)/CrsMatrices/Galeri_Pei.h               \
	 $(srcdir)/CrsMatrices/Galeri_JordanBlock.h       \
	 $(srcdir)/CrsMatrices/Galeri_Recirc2D.h          \
	 $(srcdir)/CrsMatrices/Galeri_BentPipe2D.h        \
	 $(srcdir)/FiniteElements/Galeri_AbstractGrid.h           \
	 $(srcdir)/FiniteElements/Galeri_AbstractProblem.h        \
	 $(srcdir)/FiniteElements/Galeri_AbstractQuadrature.h     \
	 $(srcdir)/FiniteElements/Galeri_AbstractVariational.h    \
	 $(srcdir)/FiniteElements/Galeri_FileGrid.h               \
	 $(srcdir)/FiniteElements/Galeri_TriangleRectangleGrid.h  \
	 $(srcdir)/FiniteElements/Galeri_QuadRectangleGrid.h      \
	 $(srcdir)/FiniteElements/Galeri_HexCubeGrid.h            \
	 $(srcdir)/FiniteElements/Galeri_TetCubeGrid.h            \
	 $(srcdir)/FiniteElements/Galeri_TRIANGLEGrid.h           \
	 $(srcdir)/FiniteElements/Galeri_QuadQuadrature.h         \
	 $(srcdir)/FiniteElements/Galeri_TetQuadrature.h          \
	 $(srcdir)/FiniteElements/Galeri_TriangleQuadrature.h     \
	 $(srcdir)/FiniteElements/Galeri_HexQuadrature.h          \
	 $(srcdir)/FiniteElements/Galeri_LinearProblem.h          \
	 $(srcdir)/FiniteElements/Galeri_GalerkinVariational.h    \
	 $(srcdir)/FiniteElements/Galeri_SUPGVariational.h        \
	 $(srcdir)/FiniteElements/Galeri_MEDITInterface.h         \
	 $(srcdir)/FiniteElements/Galeri_Workspace.h              


CORE = $(srcdir)/Galeri_Maps.cpp        \
       $(srcdir)/Galeri_CrsMatrices.cpp \
       $(srcdir)/Galeri_VbrMatrices.cpp \
       $(srcdir)/Galeri_Utils.cpp       \
       $(srcdir)/Galeri_ReadHB.cpp      \
       $(srcdir)/iohb.cpp

AM_CPPFLAGS = -I$(srcdir)                      \
              -I$(top_srcdir)/../epetra/src    \
	      -I$(top_builddir)/../epetra/src  \
	      -I$(top_srcdir)/../teuchos/src   \
	      -I$(top_builddir)/../teuchos/src

# ------------------------------------------------------------------------
# Galeri library specifications
# ------------------------------------------------------------------------
GALERI_LIB = libgaleri.a

GALERI_H = $(CORE_H) 

libgaleri_a_SOURCES = $(CORE) 

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR
libgaleri_a_AR = $(ALTERNATE_AR)
else
libgaleri_a_AR = $(AR) cru
endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libgaleri_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(GALERI_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(GALERI_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = Galeri_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Library Target for Circular Dependencies
# ------------------------------------------------------------------------

all-libs:
	cd $(top_builddir)/../epetra ; $(MAKE) all-libs
	cd $(top_builddir)/../teuchos ; $(MAKE) all-libs
	$(MAKE) libgaleri.a
