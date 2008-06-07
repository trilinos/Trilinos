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

all-local: FourByFour.xml make-links

#
# Actual programs
#

TEST_SINGLE_STRATIMIKOS_SOLVER_DRIVER = test_single_stratimikos_solver_driver
#TEST_STRATIMIKOS_DRIVER = test_stratimikos_driver

#noinst_PROGRAMS = $(TEST_SINGLE_STRATIMIKOS_SOLVER_DRIVER) $(TEST_STRATIMIKOS_DRIVER)
noinst_PROGRAMS = $(TEST_SINGLE_STRATIMIKOS_SOLVER_DRIVER)

test_single_stratimikos_solver_driver_SOURCES = $(srcdir)/test_single_stratimikos_solver_driver.cpp $(common_sources)
test_single_stratimikos_solver_driver_DEPENDENCIES = $(common_dependencies)
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
	cp  $(srcdir)/_FourByFour.aztecoo.ifpack.adjoint.xml FourByFour.aztecoo.ifpack.adjoint.xml
	cp  $(srcdir)/_FourByFour.aztecoo.ml.xml FourByFour.aztecoo.ml.xml
	cp  $(srcdir)/_FourByFour.belos.np.xml FourByFour.belos.np.xml
	cp  $(srcdir)/_FourByFour.belos.ifpack.xml FourByFour.belos.ifpack.xml
	cp  $(srcdir)/_FourByFour.belos.ml.xml FourByFour.belos.ml.xml

#
# Links
#

# Links for amesos executables

if ENABLE_STRATIMIKOS_AMESOS
LINK_TESTDRIVER_AMESOS=test_single_stratimikos_solver_driver_amesos.exe
endif

test_single_stratimikos_solver_driver_amesos.exe: test_single_stratimikos_solver_driver.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

# Links for aztecoo executables

if ENABLE_STRATIMIKOS_AZTECOO
if ENABLE_STRATIMIKOS_IFPACK
LINK_TESTDRIVER_AZTECOO_IFPACK=test_single_stratimikos_solver_driver_aztecoo_ifpack.exe
endif
if ENABLE_STRATIMIKOS_ML
LINK_TESTDRIVER_AZTECOO_ML=test_single_stratimikos_solver_driver_aztecoo_ml.exe
endif
endif

test_single_stratimikos_solver_driver_aztecoo_ifpack.exe: test_single_stratimikos_solver_driver.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

test_single_stratimikos_solver_driver_aztecoo_ml.exe: test_single_stratimikos_solver_driver.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

# Links for belos executables

if ENABLE_STRATIMIKOS_BELOS
if ENABLE_STRATIMIKOS_IFPACK
LINK_TESTDRIVER_BELOS_IFPACK=test_single_stratimikos_solver_driver_belos_ifpack.exe
endif
if ENABLE_STRATIMIKOS_ML
LINK_TESTDRIVER_BELOS_ML=test_single_stratimikos_solver_driver_belos_ml.exe
endif
endif

test_single_stratimikos_solver_driver_belos_ifpack.exe: test_single_stratimikos_solver_driver.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

test_single_stratimikos_solver_driver_belos_ml.exe: test_single_stratimikos_solver_driver.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

# Main links

make-links: \
  $(LINK_TESTDRIVER_AMESOS) \
  $(LINK_TESTDRIVER_AZTECOO_IFPACK) \
  $(LINK_TESTDRIVER_AZTECOO_ML) \
  $(LINK_TESTDRIVER_BELOS_IFPACK) \
  $(LINK_TESTDRIVER_BELOS_ML)

#
# Common stuff
#

common_dependencies = $(top_builddir)/src/libstratimikos.a

common_sources = $(srcdir)/test_single_stratimikos_solver.cpp

include $(top_builddir)/Makefile.export.stratimikos

AM_CPPFLAGS =  -I$(top_srcdir)/../thyra/example/operator_solve $(STRATIMIKOS_INCLUDES)

common_ldadd = $(STRATIMIKOS_LIBS)

EXTRA_DIST = $(srcdir)/test_single_stratimikos_solver.hpp \
	$(srcdir)/_FourByFour.xml \
	$(srcdir)/_FourByFour.amesos.xml \
	$(srcdir)/_FourByFour.aztecoo.ifpack.xml \
        $(srcdir)/_FourByFour.aztecoo.ifpack.adjoint.xml \
	$(srcdir)/_FourByFour.aztecoo.ml.xml \
	$(srcdir)/_FourByFour.aztecoo.np.xml \
	$(srcdir)/_FourByFour.belos.ifpack.xml \
	$(srcdir)/_FourByFour.belos.ml.xml \
	$(srcdir)/_FourByFour.belos.np.xml \
	$(srcdir)/definition

# This is already added part of STRATIMIKOS_LIBS and therefore automake does not need to add this again!
LIBS =
