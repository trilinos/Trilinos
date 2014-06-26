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
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA
# Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
#
# ***********************************************************************
#@HEADER

EXEEXT = .exe


bin_PROGRAMS = Rythmos_BackwardEuler_ConvergenceTest \
							 Rythmos_ForwardEuler_ConvergenceTest \
							 Rythmos_ImplicitBDF_ConvergenceTest \
							 Rythmos_ExplicitRK_ConvergenceTest \
							 Rythmos_ImplicitRK_ConvergenceTest 

#
# Actual programs
#


# ------ Backward Euler ------
Rythmos_BackwardEuler_ConvergenceTest_INCLUDES =\
  $(srcdir)/Rythmos_ConvergenceTestHelpers.hpp\
	$(srcdir)/../UnitTest/Rythmos_UnitTestModels.hpp\
  $(srcdir)/Rythmos_BackwardEuler_ConvergenceTest.hpp
Rythmos_BackwardEuler_ConvergenceTest_SOURCES =\
  $(top_srcdir)/../epetraext/example/model_evaluator/DiagonalTransient/EpetraExt_DiagonalTransientModel.cpp\
  $(srcdir)/../SinCos/SinCosModel.cpp\
	$(srcdir)/Rythmos_ConvergenceTest.cpp\
  $(srcdir)/Rythmos_ConvergenceTestHelpers.cpp\
	$(srcdir)/Rythmos_BackwardEuler_ConvergenceTest.cpp
Rythmos_BackwardEuler_ConvergenceTest_DEPENDENCIES = $(common_dependencies)
Rythmos_BackwardEuler_ConvergenceTest_LDADD = $(common_ldadd)


# ------ Forward Euler ------
Rythmos_ForwardEuler_ConvergenceTest_INCLUDES =\
  $(srcdir)/Rythmos_ConvergenceTestHelpers.hpp\
	$(srcdir)/../UnitTest/Rythmos_UnitTestModels.hpp\
  $(srcdir)/Rythmos_ForwardEuler_ConvergenceTest.hpp
Rythmos_ForwardEuler_ConvergenceTest_SOURCES =\
  $(top_srcdir)/../epetraext/example/model_evaluator/DiagonalTransient/EpetraExt_DiagonalTransientModel.cpp\
  $(srcdir)/../SinCos/SinCosModel.cpp\
	$(srcdir)/Rythmos_ConvergenceTest.cpp\
  $(srcdir)/Rythmos_ConvergenceTestHelpers.cpp\
	$(srcdir)/Rythmos_ForwardEuler_ConvergenceTest.cpp
Rythmos_ForwardEuler_ConvergenceTest_DEPENDENCIES = $(common_dependencies)
Rythmos_ForwardEuler_ConvergenceTest_LDADD = $(common_ldadd)


# ------ Implicit BDF ------
Rythmos_ImplicitBDF_ConvergenceTest_INCLUDES =\
  $(srcdir)/Rythmos_ConvergenceTestHelpers.hpp\
	$(srcdir)/../UnitTest/Rythmos_UnitTestModels.hpp\
  $(srcdir)/Rythmos_ImplicitBDF_ConvergenceTest.hpp
Rythmos_ImplicitBDF_ConvergenceTest_SOURCES =\
  $(top_srcdir)/../epetraext/example/model_evaluator/DiagonalTransient/EpetraExt_DiagonalTransientModel.cpp\
  $(srcdir)/../SinCos/SinCosModel.cpp\
	$(srcdir)/Rythmos_ConvergenceTest.cpp\
  $(srcdir)/Rythmos_ConvergenceTestHelpers.cpp\
	$(srcdir)/Rythmos_ImplicitBDF_ConvergenceTest.cpp
Rythmos_ImplicitBDF_ConvergenceTest_DEPENDENCIES = $(common_dependencies)
Rythmos_ImplicitBDF_ConvergenceTest_LDADD = $(common_ldadd)


# ------ Explicit RK ------
Rythmos_ExplicitRK_ConvergenceTest_INCLUDES =\
  $(srcdir)/Rythmos_ConvergenceTestHelpers.hpp\
	$(srcdir)/../UnitTest/Rythmos_UnitTestModels.hpp\
  $(srcdir)/Rythmos_ExplicitRK_ConvergenceTest.hpp
Rythmos_ExplicitRK_ConvergenceTest_SOURCES =\
  $(top_srcdir)/../epetraext/example/model_evaluator/DiagonalTransient/EpetraExt_DiagonalTransientModel.cpp\
  $(srcdir)/../SinCos/SinCosModel.cpp\
	$(srcdir)/Rythmos_ConvergenceTest.cpp\
  $(srcdir)/Rythmos_ConvergenceTestHelpers.cpp\
	$(srcdir)/Rythmos_ExplicitRK_ConvergenceTest.cpp
Rythmos_ExplicitRK_ConvergenceTest_DEPENDENCIES = $(common_dependencies)
Rythmos_ExplicitRK_ConvergenceTest_LDADD = $(common_ldadd)


# ------ Implicit RK ------
Rythmos_ImplicitRK_ConvergenceTest_INCLUDES =\
  $(srcdir)/Rythmos_ConvergenceTestHelpers.hpp\
	$(srcdir)/../UnitTest/Rythmos_UnitTestModels.hpp\
  $(srcdir)/Rythmos_ImplicitRK_ConvergenceTest.hpp
Rythmos_ImplicitRK_ConvergenceTest_SOURCES =\
  $(top_srcdir)/../epetraext/example/model_evaluator/DiagonalTransient/EpetraExt_DiagonalTransientModel.cpp\
  $(srcdir)/../SinCos/SinCosModel.cpp\
	$(srcdir)/Rythmos_ConvergenceTest.cpp\
  $(srcdir)/Rythmos_ConvergenceTestHelpers.cpp\
	$(srcdir)/Rythmos_ImplicitRK_ConvergenceTest.cpp
Rythmos_ImplicitRK_ConvergenceTest_DEPENDENCIES = $(common_dependencies)
Rythmos_ImplicitRK_ConvergenceTest_LDADD = $(common_ldadd)


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
