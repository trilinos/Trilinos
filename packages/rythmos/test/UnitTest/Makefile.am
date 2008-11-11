#@HEADER
# ***********************************************************************
#
#                           Rythmos Package
#                 Copyright (2006) Sandia Corporation
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
# Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
#
# ***********************************************************************
#@HEADER

EXEEXT = .exe


bin_PROGRAMS = Rythmos_UnitTest

#
# Actual programs
#

EXTRA_DIST =\
  $(srcdir)/Rythmos_UnitTestHelpers.hpp \
	$(srcdir)/Rythmos_UnitTestModels.hpp \
	$(srcdir)/../SinCos/SinCosModel.hpp \
	$(srcdir)/../PolynomialModel/PolynomialModel.hpp

Rythmos_UnitTest_SOURCES =\
  $(top_srcdir)/../epetraext/example/model_evaluator/DiagonalTransient/EpetraExt_DiagonalTransientModel.cpp\
	$(srcdir)/Rythmos_ConvergenceTestHelpers_UnitTest.cpp\
  $(srcdir)/Rythmos_DataStore_UnitTest.cpp\
  $(srcdir)/Rythmos_ExplicitRK_UnitTest.cpp\
  $(srcdir)/Rythmos_ImplicitBDF_UnitTest.cpp\
  $(srcdir)/Rythmos_ImplicitRK_UnitTest.cpp\
  $(srcdir)/Rythmos_InterpolationBuffer_UnitTest.cpp\
	$(srcdir)/Rythmos_LinearInterpolator_UnitTest.cpp\
	$(srcdir)/Rythmos_Quadrature_UnitTest.cpp\
  $(srcdir)/Rythmos_RKButcherTableau_UnitTest.cpp\
  $(srcdir)/../SinCos/SinCosModel.cpp\
  $(srcdir)/../PolynomialModel/PolynomialModel.cpp\
  $(srcdir)/Rythmos_SinCosModel_UnitTest.cpp\
  $(srcdir)/Rythmos_TimeRange_UnitTest.cpp\
	$(srcdir)/Rythmos_UnitTest.cpp\
	$(srcdir)/Rythmos_UnitTestHelpers.cpp
Rythmos_UnitTest_DEPENDENCIES = $(common_dependencies)
Rythmos_UnitTest_LDADD = $(common_ldadd)

#
# Common stuff
# 

common_dependencies = \
  $(top_builddir)/src/librythmos.a

include $(top_builddir)/Makefile.export.rythmos

_EXPORT_INCLUDES = $(RYTHMOS_TEST_INCLUDES)
_EXPORT_LIBS = $(RYTHMOS_TEST_LIBS)

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(_EXPORT_INCLUDES))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(_EXPORT_LIBS))
else
EXPORT_INCLUDES = $(_EXPORT_INCLUDES)
EXPORT_LIBS = $(_EXPORT_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

common_ldadd = $(EXPORT_LIBS)

CPPFLAGS += -I$(top_srcdir)/../epetraext/example/model_evaluator/DiagonalTransient

# This is already added part of RYTHMOS_LIBS and therefore automake does not need to add this again!
LIBS =
