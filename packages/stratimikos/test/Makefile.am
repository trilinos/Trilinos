# @HEADER
# ************************************************************************
# 
#         Stratimikos: Thyra-based strategies for linear solvers
#                Copyright (2006) Sandia Corporation
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
# Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
# 
# ************************************************************************
# @HEADER

EXEEXT = .exe

#
# Actual programs
#

TEST_SINGLE_STRATIMIKOS_SOLVER_DRIVER = test_single_stratimikos_solver_driver
#TEST_STRATIMIKOS_DRIVER = test_stratimikos_driver

#noinst_PROGRAMS = $(TEST_SINGLE_STRATIMIKOS_SOLVER_DRIVER) $(TEST_STRATIMIKOS_DRIVER)
noinst_PROGRAMS = $(TEST_SINGLE_STRATIMIKOS_SOLVER_DRIVER)

test_single_stratimikos_solver_driver_SOURCES = $(srcdir)/test_single_stratimikos_solver_driver.cpp $(common_sources)
test_single_stratimikos_solver_driver_DEPENDENCIES = $(common_dependencies)  FourByFour.xml
test_single_stratimikos_solver_driver_LDADD = $(common_ldadd)

#test_stratimikos_driver_SOURCES = $(srcdir)/test_stratimikos_driver.cpp $(common_sources)
#test_stratimikos_driver_DEPENDENCIES = $(common_dependencies) FourByFour.mtx
#test_stratimikos_driver_LDADD = $(common_ldadd)

FourByFour.xml:
	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_FourByFour.mtx FourByFour.mtx                          
	cp  $(srcdir)/_FourByFour.xml FourByFour.xml
	cp  $(srcdir)/_FourByFour.amesos.xml FourByFour.amesos.xml
	cp  $(srcdir)/_FourByFour.aztecoo.np.xml FourByFour.aztecoo.np.xml
	cp  $(srcdir)/_FourByFour.aztecoo.ifpack.xml FourByFour.aztecoo.ifpack.xml
	cp  $(srcdir)/_FourByFour.belos.np.xml FourByFour.belos.np.xml
	cp  $(srcdir)/_FourByFour.belos.ifpack.xml FourByFour.belos.ifpack.xml
#	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_bcsstk01.mtx bcsstk01.mtx                              
#	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_bcsstk02.mtx bcsstk02.mtx                              
#	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_bcsstk04.mtx bcsstk04.mtx                              
#	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_nos1.mtx nos1.mtx                                      
#	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_nos5.mtx nos5.mtx                                      
#	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_Diagonal.mtx Diagonal.mtx                              
#	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_MissingADiagonal.mtx MissingADiagonal.mtx                              
#	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_KheadK.mtx KheadK.mtx                                  
#	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_KheadSorted.mtx  KheadSorted.mtx
#	cp  $(top_srcdir)/../amesos/test/Test_Basic/In_Khead.triS Khead.triS  


#
# Common stuff
#

common_sources = $(srcdir)/test_single_stratimikos_solver.cpp

common_dependencies = $(top_builddir)/src/libstratimikos.a

include $(top_builddir)/Makefile.export.stratimikos

AM_CPPFLAGS =  -I$(top_srcdir)/../thyra/example/operator_solve $(STRATIMIKOS_INCLUDES)

common_ldadd = $(STRATIMIKOS_LIBS)

EXTRA_DIST = $(srcdir)/test_single_stratimikos_solver.hpp

# This is already added part of STRATIMIKOS_LIBS and therefore automake does not need to add this again!
LIBS =
