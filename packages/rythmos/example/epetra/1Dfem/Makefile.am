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


#
# Common stuff
#

include $(top_builddir)/Makefile.export.rythmos
include $(top_builddir)/../nox/Makefile.export.nox

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

EXEEXT = .exe

all-local: $(1DFEM_EXE_LINKS)
	@echo ""
	@echo "rythmos/example/epetra/1Dfem/1DfemTransient built successfully."
	@echo ""

#
# Actual programs
#

noinst_PROGRAMS = 1DfemTransient

1DfemTransient_SOURCES = \
	                 $(srcdir)/1DfemInterface.H \
                   $(srcdir)/1DfemInterface.C \
                   $(srcdir)/1DfemTransient.H \
									 $(srcdir)/ExampleApplication1Dfem.hpp \
									 $(srcdir)/ExampleApplication1Dfem.cpp \
									 $(srcdir)/cxx_main.cpp

1DfemTransient_DEPENDENCIES = \
                              $(top_builddir)/../nox/src/libnox.a \
                              $(top_builddir)/../nox/src-epetra/libnoxepetra.a \
                              $(top_builddir)/src/librythmos.a

1DfemTransient_LDADD = $(EXPORT_LIBS)


#
# Links
#

if ENABLE_RYTHMOS_STRATIMIKOS_AMESOS
LINK_DIAGONAL_AMESOS=1DfemTransient_amesos.exe
endif

if ENABLE_RYTHMOS_STRATIMIKOS_AZTECOO
LINK_DIAGONAL_AZTECOO=1DfemTransient_aztecoo.exe
if ENABLE_RYTHMOS_STRATIMIKOS_IFPACK
LINK_DIAGONAL_AZTECOO_IFPACK=1DfemTransient_aztecoo_ifpack.exe
endif
if ENABLE_RYTHMOS_STRATIMIKOS_ML
LINK_DIAGONAL_AZTECOO_ML=1DfemTransient_aztecoo_ml.exe
endif
endif

if ENABLE_RYTHMOS_STRATIMIKOS_BELOS
LINK_DIAGONAL_BELOS=1DfemTransient_belos.exe
if ENABLE_RYTHMOS_STRATIMIKOS_IFPACK
LINK_DIAGONAL_BELOS_IFPACK=1DfemTransient_belos_ifpack.exe
endif
if ENABLE_RYTHMOS_STRATIMIKOS_ML
LINK_DIAGONAL_BELOS_ML=1DfemTransient_belos_ml.exe
endif
endif

1DfemTransient_amesos.exe: 1DfemTransient.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

1DfemTransient_aztecoo.exe: 1DfemTransient.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

1DfemTransient_aztecoo_ifpack.exe: 1DfemTransient.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

1DfemTransient_aztecoo_ml.exe: 1DfemTransient.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

1DfemTransient_belos.exe: 1DfemTransient.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

1DfemTransient_belos_ifpack.exe: 1DfemTransient.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

1DfemTransient_belos_ml.exe: 1DfemTransient.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

1DFEM_EXE_LINKS = \
  $(LINK_DIAGONAL_AMESOS) \
  $(LINK_DIAGONAL_AZTECOO) \
  $(LINK_DIAGONAL_AZTECOO_IFPACK) \
  $(LINK_DIAGONAL_AZTECOO_ML) \
  $(LINK_DIAGONAL_BELOS) \
  $(LINK_DIAGONAL_BELOS_IFPACK) \
  $(LINK_DIAGONAL_BELOS_ML)
