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

all-local: generate-sample-output.pl amesos.klu.xml MixedOrderPhysicsBasedPreconditioner.Base.xml make-links

generate-sample-output.pl : generate-sample-output.pl.stub
	cp generate-sample-output.pl.stub generate-sample-output.pl
	chmod a+x generate-sample-output.pl

#
# Actual programs
#

SIMPLE_STRATIMIKOS_EXAMPLE = simple_stratimikos_example
MIXED_ORDER_PHYSICS_BASED_PRECONDITIONER = MixedOrderPhysicsBasedPreconditioner

noinst_PROGRAMS = $(SIMPLE_STRATIMIKOS_EXAMPLE) $(MIXED_ORDER_PHYSICS_BASED_PRECONDITIONER)

simple_stratimikos_example_SOURCES = $(srcdir)/simple_stratimikos_example.cpp
simple_stratimikos_example_DEPENDENCIES = $(common_dependencies)
simple_stratimikos_example_LDADD = $(common_ldadd)

amesos.klu.xml:
	cp $(srcdir)/_amesos.klu.xml amesos.klu.xml
	cp $(srcdir)/_aztecoo.ifpack.xml aztecoo.ifpack.xml
	cp $(srcdir)/_aztecoo.ifpack.xml aztecoo.ml.xml
	cp $(srcdir)/_belos.ifpack.xml belos.ifpack.xml
	cp $(srcdir)/_belos.ifpack.xml belos.ml.xml
	cp $(top_srcdir)/../amesos/test/Test_Basic/In_FourByFour.mtx FourByFour.mtx

MixedOrderPhysicsBasedPreconditioner_SOURCES = $(srcdir)/MixedOrderPhysicsBasedPreconditioner.cpp
MixedOrderPhysicsBasedPreconditioner_DEPENDENCIES = $(common_dependencies)
MixedOrderPhysicsBasedPreconditioner_LDADD = $(common_ldadd)

MixedOrderPhysicsBasedPreconditioner.Base.xml:
	cp $(srcdir)/_MixedOrderPhysicsBasedPreconditioner.Base.xml MixedOrderPhysicsBasedPreconditioner.Base.xml
	cp $(srcdir)/_MixedOrderPhysicsBasedPreconditioner.AllSolveOutput.xml MixedOrderPhysicsBasedPreconditioner.AllSolveOutput.xml
	cp $(srcdir)/_MixedOrderPhysicsBasedPreconditioner.P2SolveOutputMSolveOutputToFile.xml MixedOrderPhysicsBasedPreconditioner.P2SolveOutputMSolveOutputToFile.xml
	cp $(srcdir)/_MixedOrderPhysicsBasedPreconditioner.P2SolveOutput.xml MixedOrderPhysicsBasedPreconditioner.P2SolveOutput.xml
	cp $(srcdir)/_MixedOrderPhysicsBasedPreconditioner.Belos.xml MixedOrderPhysicsBasedPreconditioner.Belos.xml
	cp $(srcdir)/_P2.mtx P2.mtx
	cp $(srcdir)/_P1.mtx P1.mtx
	cp $(srcdir)/_M11.mtx M11.mtx
	cp $(srcdir)/_M12.mtx M12.mtx
	cp $(srcdir)/_M21.mtx M21.mtx
	cp $(srcdir)/_M22.mtx M22.mtx


#
# Links
#


# Links for amesos executables

if ENABLE_STRATIMIKOS_AMESOS
LINK_SIMPLE_AMESOS=simple_stratimikos_example_amesos.exe
LINK_MIXEDORDER_AMESOS=MixedOrderPhysicsBasedPreconditioner_amesos.exe
endif

simple_stratimikos_example_amesos.exe: simple_stratimikos_example.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

MixedOrderPhysicsBasedPreconditioner_amesos.exe: MixedOrderPhysicsBasedPreconditioner.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

# Links for aztecoo executables

if ENABLE_STRATIMIKOS_AZTECOO
if ENABLE_STRATIMIKOS_IFPACK
LINK_SIMPLE_AZTECOO_IFPACK=simple_stratimikos_example_aztecoo_ifpack.exe
LINK_MIXEDORDER_AZTECOO_IFPACK=MixedOrderPhysicsBasedPreconditioner_aztecoo_ifpack.exe
endif
if ENABLE_STRATIMIKOS_ML
LINK_SIMPLE_AZTECOO_ML=simple_stratimikos_example_aztecoo_ml.exe
LINK_MIXEDORDER_AZTECOO_ML=MixedOrderPhysicsBasedPreconditioner_aztecoo_ml.exe
endif
endif

simple_stratimikos_example_aztecoo_ifpack.exe: simple_stratimikos_example.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

MixedOrderPhysicsBasedPreconditioner_aztecoo_ifpack.exe: MixedOrderPhysicsBasedPreconditioner.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

simple_stratimikos_example_aztecoo_ml.exe: simple_stratimikos_example.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

MixedOrderPhysicsBasedPreconditioner_aztecoo_ml.exe: MixedOrderPhysicsBasedPreconditioner.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

# Links for belos executables

if ENABLE_STRATIMIKOS_BELOS
if ENABLE_STRATIMIKOS_IFPACK
LINK_SIMPLE_BELOS_IFPACK=simple_stratimikos_example_belos_ifpack.exe
LINK_MIXEDORDER_BELOS_IFPACK=MixedOrderPhysicsBasedPreconditioner_belos_ifpack.exe
endif
if ENABLE_STRATIMIKOS_ML
LINK_SIMPLE_BELOS_ML=simple_stratimikos_example_belos_ml.exe
LINK_MIXEDORDER_BELOS_ML=MixedOrderPhysicsBasedPreconditioner_belos_ml.exe
endif
endif

simple_stratimikos_example_belos_ifpack.exe: simple_stratimikos_example.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

MixedOrderPhysicsBasedPreconditioner_belos_ifpack.exe: MixedOrderPhysicsBasedPreconditioner.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

simple_stratimikos_example_belos_ml.exe: simple_stratimikos_example.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

MixedOrderPhysicsBasedPreconditioner_belos_ml.exe: MixedOrderPhysicsBasedPreconditioner.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

# Main links

make-links: \
  $(LINK_SIMPLE_AMESOS) $(LINK_MIXEDORDER_AMESOS) \
  $(LINK_SIMPLE_AZTECOO_IFPACK) $(LINK_MIXEDORDER_AZTECOO_IFPACK) \
  $(LINK_SIMPLE_AZTECOO_ML) $(LINK_MIXEDORDER_AZTECOO_ML) \
  $(LINK_SIMPLE_BELOS_IFPACK) $(LINK_MIXEDORDER_BELOS_IFPACK) \
  $(LINK_SIMPLE_BELOS_ML) $(LINK_MIXEDORDER_BELOS_ML)

#
# Common stuff
#

common_dependencies = $(top_builddir)/src/libstratimikos.a

include $(top_builddir)/Makefile.export.stratimikos

AM_CPPFLAGS =  $(STRATIMIKOS_INCLUDES)

common_ldadd = $(STRATIMIKOS_LIBS)

EXTRA_DIST = $(srcdir)/simple_stratimikos_example.cpp \
	$(srcdir)/_amesos.klu.xml \
	$(srcdir)/_aztecoo.ifpack.xml \
	$(srcdir)/_belos.ifpack.xml \
	$(srcdir)/MixedOrderPhysicsBasedPreconditioner.cpp \
	$(srcdir)/_MixedOrderPhysicsBasedPreconditioner.Base.xml \
	$(srcdir)/_MixedOrderPhysicsBasedPreconditioner.AllSolveOutput.xml \
	$(srcdir)/_MixedOrderPhysicsBasedPreconditioner.P2SolveOutputMSolveOutputToFile.xml \
	$(srcdir)/_MixedOrderPhysicsBasedPreconditioner.P2SolveOutput.xml \
	$(srcdir)/_MixedOrderPhysicsBasedPreconditioner.Belos.xml \
	$(srcdir)/_P2.mtx \
	$(srcdir)/_P1.mtx \
	$(srcdir)/_M11.mtx \
	$(srcdir)/_M12.mtx \
	$(srcdir)/_M21.mtx \
	$(srcdir)/_M22.mtx

# This is already added part of STRATIMIKOS_LIBS and therefore automake does not need to add this again!
LIBS =
