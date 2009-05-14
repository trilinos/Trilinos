## $Id$
## $Source$
# @HEADER
# ************************************************************************
# 
#                     Stokhos Package
#                 Copyright (2008) Sandia Corporation
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
# Questions? Contact Eric T. Phipps (etphipp@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package stokhos subdirectory src built successfully."
	@echo ""

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside Stokhos. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

STOKHOS_H = \
	$(srcdir)/Stokhos_ConfigDefs.h \
	$(srcdir)/Stokhos_Version.hpp \
	$(srcdir)/Stokhos_gaussq.h \
	$(srcdir)/Stokhos_DynamicArrayTraits.hpp \
	$(srcdir)/Stokhos_Polynomial.hpp \
	$(srcdir)/Stokhos_PolynomialImp.hpp \
	$(srcdir)/Stokhos_OneDOrthogPolyBasis.hpp \
	$(srcdir)/Stokhos_OneDOrthogPolyBasisBase.hpp \
	$(srcdir)/Stokhos_OneDOrthogPolyBasisBaseImp.hpp \
	$(srcdir)/Stokhos_HermiteBasis.hpp \
	$(srcdir)/Stokhos_HermiteBasisImp.hpp \
	$(srcdir)/Stokhos_LegendreBasis.hpp \
	$(srcdir)/Stokhos_LegendreBasisImp.hpp \
	$(srcdir)/Stokhos_OrthogPolyBasis.hpp \
	$(srcdir)/Stokhos_CompletePolynomialBasis.hpp \
	$(srcdir)/Stokhos_CompletePolynomialBasisImp.hpp \
	$(srcdir)/Stokhos_Sparse3Tensor.hpp \
	$(srcdir)/Stokhos_Sparse3TensorImp.hpp \
	$(srcdir)/Stokhos_Dense3Tensor.hpp \
	$(srcdir)/Stokhos_Dense3TensorImp.hpp \
	$(srcdir)/Stokhos_Quadrature.hpp \
	$(srcdir)/Stokhos_TensorProductQuadrature.hpp \
	$(srcdir)/Stokhos_TensorProductQuadratureImp.hpp \
	$(srcdir)/Stokhos_SparseGridQuadrature.hpp \
	$(srcdir)/Stokhos_SparseGridQuadratureImp.hpp \
	$(srcdir)/Stokhos_OrthogPolyApprox.hpp \
	$(srcdir)/Stokhos_OrthogPolyApproxImp.hpp \
	$(srcdir)/Stokhos_OrthogPolyExpansion.hpp \
	$(srcdir)/Stokhos_QuadOrthogPolyExpansion.hpp \
	$(srcdir)/Stokhos_QuadOrthogPolyExpansionImp.hpp \
	$(srcdir)/Stokhos_DerivOrthogPolyExpansion.hpp \
	$(srcdir)/Stokhos_DerivOrthogPolyExpansionImp.hpp \
	$(srcdir)/Stokhos_ForUQTKOrthogPolyExpansion.hpp \
	$(srcdir)/Stokhos_ForUQTKOrthogPolyExpansionImp.hpp \
	$(srcdir)/Stokhos_MatrixFreeEpetraOp.hpp \
	$(srcdir)/Stokhos_MeanEpetraOp.hpp

STOKHOS = \
	$(srcdir)/Stokhos_Version.cpp \
	$(srcdir)/Stokhos_MatrixFreeEpetraOp.cpp \
	$(srcdir)/Stokhos_MeanEpetraOp.cpp

# Add packages that this package is dependent on
include $(top_builddir)/Makefile.export.stokhos

# ------------------------------------------------------------------------
# Stokhos library specifications
# ------------------------------------------------------------------------

AM_CPPFLAGS = $(STOKHOS_INCLUDES)

STOKHOS_LIB = libstokhos.a

libstokhos_a_SOURCES = $(STOKHOS) $(STOKHOS_H)

EXTRA_libstokhos_a_SOURCES = 

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

#np# replace stokhos with the name of the package being autotool'ed here
libstokhos_a_AR = $(ALTERNATE_AR)
else

#np# replace stokhos with the name of the package being autotool'ed here
libstokhos_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libstokhos_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(STOKHOS_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(STOKHOS_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Stokhos_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include

nodist_execinclude_HEADERS = Stokhos_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






