#@HEADER
# ************************************************************************
# 
#                    Teuchos: Common Tools Package
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
#@HEADER

EXEEXT = .exe

include $(top_builddir)/Makefile.export.teuchos

if HAVE_GCC_CHECKEDSTL
  GCC_CHECKEDSTL_EXE_NAME=GCC_CheckedSTL_UnitTests
endif


bin_PROGRAMS = RCP_test Ptr_test ArrayRCP_test Array_test \
  ArrayView_test Tuple_test MemoryManagement_UnitTests $(GCC_CHECKEDSTL_EXE_NAME)

RCP_test_SOURCES = $(srcdir)/RCP_test.cpp $(common_sources)
RCP_test_DEPENDENCIES = $(top_builddir)/src/libteuchos.a
RCP_test_LDADD = $(TEUCHOS_LIBS)

Ptr_test_SOURCES = $(srcdir)/Ptr_test.cpp $(common_sources)
Ptr_test_DEPENDENCIES = $(top_builddir)/src/libteuchos.a
Ptr_test_LDADD = $(TEUCHOS_LIBS)

ArrayRCP_test_SOURCES = $(srcdir)/ArrayRCP_test.cpp
ArrayRCP_test_DEPENDENCIES = $(top_builddir)/src/libteuchos.a
ArrayRCP_test_LDADD = $(TEUCHOS_LIBS)

Array_test_SOURCES = $(srcdir)/Array_test.cpp
Array_test_DEPENDENCIES = $(top_builddir)/src/libteuchos.a
Array_test_LDADD = $(TEUCHOS_LIBS)

ArrayView_test_SOURCES = $(srcdir)/ArrayView_test.cpp
ArrayView_test_DEPENDENCIES = $(top_builddir)/src/libteuchos.a
ArrayView_test_LDADD = $(TEUCHOS_LIBS)

Tuple_test_SOURCES = $(srcdir)/Tuple_test.cpp
Tuple_test_DEPENDENCIES = $(top_builddir)/src/libteuchos.a
Tuple_test_LDADD = $(TEUCHOS_LIBS)

MemoryManagement_UnitTests_SOURCES = \
  $(srcdir)/RCP_UnitTests.cpp \
  $(srcdir)/RCPNodeHandle_UnitTests.cpp \
  $(srcdir)/Array_UnitTest_helpers.cpp \
  $(srcdir)/Array_UnitTests.cpp \
  $(srcdir)/ArrayView_UnitTests.cpp \
  $(srcdir)/ArrayConversions_UnitTests.cpp \
  $(unit_test_common_sources) \
  $(common_sources)
MemoryManagement_UnitTests_DEPENDENCIES = $(top_builddir)/src/libteuchos.a
MemoryManagement_UnitTests_LDADD = $(TEUCHOS_LIBS)

if HAVE_GCC_CHECKEDSTL
  GCC_CheckedSTL_UnitTests_SOURCES = $(srcdir)/GCC_CheckedSTL_UnitTests.cpp \
  $(common_sources)
  GCC_CheckedSTL_UnitTests_DEPENDENCIES = $(top_builddir)/src/libteuchos.a
  GCC_CheckedSTL_UnitTests_LDADD = $(TEUCHOS_LIBS)
endif

unit_test_common_sources = \
  $(top_srcdir)/test/UnitTest/Teuchos_StandardUnitTestMain.cpp

common_sources = \
   $(srcdir)/TestClasses.cpp

AM_CPPFLAGS = $(TEUCHOS_INCLUDES)

# This is already added part of TEUCHOS_LIBS and therefore automake does not need to add this again!
LIBS =


EXTRA_DIST = TestClasses.hpp Array_UnitTest_helpers.hpp
