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

# Break tradoptest.c0 into 64 pieces for use with g++ -O3, which uses
# too much memory if we try to compile all tests at once.

# Include export Makefiles for correct paths
include $(top_builddir)/Makefile.export.sacado

EXEEXT = .exe

noinst_PROGRAMS = hesopcheck hesopcheck_ADC

testall : $(noinst_PROGRAMS) tradoptest.exe
	./hesopcheck.exe
	./hesopcheck_ADC.exe
	for i in `ls tradoptest_*.exe`; do ./$$i; done

INCLUDES = $(SACADO_INCLUDES)
hesopcheck_SOURCES = hesopcheck.cpp
hesopcheck_CXXFLAGS =
hesopcheck_DEPENDENCIES = $(top_builddir)/src/libsacado.a tradoptest.exe
hesopcheck_LDADD = $(SACADO_LIBS)
hesopcheck_ADC_SOURCES = hesopcheck.cpp $(top_builddir)/src/Sacado_radops2.cpp
hesopcheck_ADC_CXXFLAGS = -DRAD_AUTO_AD_Const
hesopcheck_ADC_DEPENDENCIES = $(top_builddir)/src/libsacado.a tradoptest.exe
hesopcheck_ADC_LDADD = $(SACADO_LIBS)

# If we are configured with --with-gnumake, use GNU's extensions to make
# to compile all of the tradoptests a lot faster if using -j xx for
# parallel make
if USING_GNUMAKE

export CXXCOMPILE

tradoptest.exe : testgen.as tradoptest.c0
	$(SHELL) testgen.as $(srcdir)/tradoptest.c0
	$(MAKE) -f $(srcdir)/GNUMakefile
	echo '#!/bin/sh' >tradoptest.exe; echo 'echo OK' >>tradoptest.exe; chmod +x tradoptest.exe

else

tradoptest.exe : testgen.as tradoptest.c0
	$(SHELL) testgen.as $(srcdir)/tradoptest.c0
	set -e;\
	for i in `ls tradoptest_*.cpp | sed 's/\.cpp//'`; do \
		echo $$i: ;\
		$(CXXCOMPILE) -o $$i.exe $$i.cpp $(LIBS);\
		$(CXXCOMPILE) -o $$i"_ADC.exe" -DRAD_AUTO_AD_Const $$i.cpp $(LIBS);\
		$(CXXCOMPILE) -o $$i"_EQA.exe" -DRAD_EQ_ALIAS $$i.cpp $(LIBS);\
		$(CXXCOMPILE) -o $$i"_RE1.exe" -DRAD_REINIT=1 $$i.cpp $(LIBS);\
		$(CXXCOMPILE) -o $$i"_RE2.exe" -DRAD_REINIT=2 $$i.cpp $(LIBS);\
		done;\
	echo '#!/bin/sh' >tradoptest.exe; echo 'echo OK' >>tradoptest.exe; chmod +x tradoptest.exe

endif

CLEANFILES = tradoptest*.cpp *.exe *.o *.obj

EXTRA_DIST = tradoptest.c0 GNUMakefile
