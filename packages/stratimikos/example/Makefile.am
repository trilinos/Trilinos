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

# The following line helps the test harness recover from build errors.

all-local: generate-sample-output.pl
	@echo ""
	@echo "Trilinos package stratimikos subdirectory example built successfully."
	@echo ""

generate-sample-output.pl : generate-sample-output.pl.stub
	cp generate-sample-output.pl.stub generate-sample-output.pl
	chmod a+x generate-sample-output.pl


#
# Actual programs
#

SIMPLE_STRATIMIKOS_EXAMPLE = simple_stratimikos_example

noinst_PROGRAMS = $(SIMPLE_STRATIMIKOS_EXAMPLE)

simple_stratimikos_example_SOURCES = $(srcdir)/simple_stratimikos_example.cpp
simple_stratimikos_example_DEPENDENCIES = $(common_dependencies) amesos.klu.xml
simple_stratimikos_example_LDADD = $(common_ldadd)

amesos.klu.xml:
	cp $(srcdir)/_amesos.klu.xml amesos.klu.xml
	cp $(srcdir)/_aztecoo.ifpack.xml aztecoo.ifpack.xml
	cp $(srcdir)/_aztecoo.ifpack.xml aztecoo.ml.xml
	cp $(srcdir)/_belos.ifpack.xml belos.ifpack.xml
	cp $(srcdir)/_belos.ifpack.xml belos.ml.xml
	cp $(top_srcdir)/../amesos/test/Test_Basic/In_FourByFour.mtx FourByFour.mtx

#
# Common stuff
#

if ENABLE_STRATIMIKOS_ML_THYRA
if ENABLE_STRATIMIKOS_ML_NOX

NOX_LIBS = $(top_builddir)/../nox/src/libnox.a

$(top_builddir)/../nox/src/libnox.a:
	cd $(top_builddir)/../nox; make all-libs

endif
endif

common_dependencies = $(top_builddir)/src/libstratimikos.a $(NOX_LIBS)

include $(top_builddir)/Makefile.export.stratimikos

AM_CPPFLAGS =  $(STRATIMIKOS_INCLUDES)

common_ldadd = $(STRATIMIKOS_LIBS)

EXTRA_DIST = $(srcdir)/simple_stratimikos_example.cpp \
	$(srcdir)/_amesos.klu.xml \
	$(srcdir)/_aztecoo.ifpack.xml \
	$(srcdir)/_belos.ifpack.xml

# This is already added part of STRATIMIKOS_LIBS and therefore automake does not need to add this again!
LIBS =
