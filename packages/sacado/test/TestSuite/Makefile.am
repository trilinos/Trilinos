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

# Include export Makefiles for correct paths
include $(top_builddir)/Makefile.export.sacado

noinst_PROGRAMS = TestSuite.exe

if ENABLE_ADOLC

TAYLOR_SRC = \
	$(srcdir)/TayUnitTests.hpp \
	$(srcdir)/TayUnitTests.cpp

TAYLOR_CXXFLAGS = -I$(ADOLC_INCDIR)

TAYLOR_LDADD = -L$(ADOLC_LIBDIR) $(ADOLC_LIBS)

endif

if ENABLE_SACADO_STOKHOS

PCE_SRC = $(srcdir)/HermiteUnitTests.hpp $(srcdir)/HermiteUnitTests.cpp

endif

if ENABLE_SACADO_TEUCHOS

BLAS_SRC = $(srcdir)/FadBLASUnitTests.hpp $(srcdir)/FadBLASUnitTests.cpp

endif

TestSuite_exe_SOURCES = \
	$(srcdir)/FadUnitTests.hpp \
	$(srcdir)/FadUnitTests2.hpp \
	$(srcdir)/FadUnitTests.cpp \
	$(srcdir)/RealFadUnitTests.cpp \
	$(srcdir)/ComplexFadUnitTests.cpp \
	$(srcdir)/ELRFadUnitTests.cpp \
	$(srcdir)/RealELRFadUnitTests.cpp \
	$(srcdir)/ComplexELRFadUnitTests.cpp \
	$(srcdir)/CacheFadUnitTests.cpp \
	$(BLAS_SRC) \
	$(TAYLOR_SRC) \
	$(PCE_SRC) \
	$(srcdir)/LogicalSparseUnitTests.hpp \
	$(srcdir)/LogicalSparseUnitTests.cpp \
	$(srcdir)/TestSuite.cpp

TestSuite_exe_CXXFLAGS = \
	$(SACADO_INCLUDES) \
	-I$(top_srcdir)/test/utils \
	-I$(top_srcdir)/Fad \
	$(TAYLOR_CXXFLAGS) \
	$(CPPUNIT_CFLAGS)

TestSuite_exe_DEPENDENCIES = \
	$(top_builddir)/test/utils/libsacadoutils.a 

TestSuite_exe_LDADD = \
	-L$(top_builddir)/test/utils -lsacadoutils \
	$(TAYLOR_LDADD) \
	$(CPPUNIT_LIBS) \
	$(SACADO_LIBS)
