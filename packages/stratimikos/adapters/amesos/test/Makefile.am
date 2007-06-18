# @HEADER
# ************************************************************************
# 
#                 Amesos: Direct Sparse Solver Package
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
# @HEADER

EXEEXT = .exe

#
# Actual programs
#

TEST_SINGLE_AMESOS_THYRA_SOLVER_DRIVER = test_single_amesos_thyra_solver_driver
TEST_AMESOS_THYRA_DRIVER = test_amesos_thyra_driver

# 2007/06/15: rabartl : ToDo: I may have to put in configure logic to only
# build these examples if epetraext is enabled!

noinst_PROGRAMS = $(TEST_SINGLE_AMESOS_THYRA_SOLVER_DRIVER) $(TEST_AMESOS_THYRA_DRIVER)

test_single_amesos_thyra_solver_driver_SOURCES = $(srcdir)/test_single_amesos_thyra_solver_driver.cpp $(common_sources)
test_single_amesos_thyra_solver_driver_DEPENDENCIES = $(common_dependencies)
test_single_amesos_thyra_solver_driver_LDADD = $(common_ldadd)

test_amesos_thyra_driver_SOURCES = $(srcdir)/test_amesos_thyra_driver.cpp $(common_sources)
test_amesos_thyra_driver_DEPENDENCIES = $(common_dependencies) FourByFour.mtx
test_amesos_thyra_driver_LDADD = $(common_ldadd)

FourByFour.mtx:
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_bcsstk01.mtx bcsstk01.mtx                              
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_bcsstk02.mtx bcsstk02.mtx                              
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_bcsstk04.mtx bcsstk04.mtx                              
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_nos1.mtx nos1.mtx                                      
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_nos5.mtx nos5.mtx                                      
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_FourByFour.mtx FourByFour.mtx                          
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_Diagonal.mtx Diagonal.mtx                              
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_MissingADiagonal.mtx MissingADiagonal.mtx                              
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_KheadK.mtx KheadK.mtx                                  
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_KheadSorted.mtx  KheadSorted.mtx
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_Khead.triS Khead.triS  

#
# Common stuff
#

common_sources = $(srcdir)/test_single_amesos_thyra_solver.cpp

common_dependencies = $(top_builddir)/adapters/amesos/src/libstratimikosamesos.a $(top_builddir)/../amesos/src/libamesos.a

include $(top_builddir)/Makefile.export.stratimikos

_EXPORT_INCLUDES = $(STRATIMIKOS_EXTRA_INCLUDES) $(STRATIMIKOS_INCLUDES)
_EXPORT_LIBS = $(STRATIMIKOS_LIBS)  $(STRATIMIKOS_EXTRA_LIBS)

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(_EXPORT_INCLUDES))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl  $(_EXPORT_LIBS))
else
EXPORT_INCLUDES = $(_EXPORT_INCLUDES)
EXPORT_LIBS = $(_EXPORT_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES) -I$(top_srcdir)/../thyra/example/operator_solve

common_ldadd = $(EXPORT_LIBS)

EXTRA_DIST = \
	$(srcdir)/test_single_amesos_thyra_solver.hpp

# This is already added part of AMESOS_LIBS and therefore automake does not
# need to add this again!
LIBS =
