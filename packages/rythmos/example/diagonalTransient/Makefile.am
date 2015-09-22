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

all-local: $(DIAGONAL_EXE_LINKS) $(XML_FILES)
	@echo ""
	@echo "rythmos/example/diagonalTransient built successfully."
	@echo ""

#
# Actual programs
#

EXAMPLE = diagonalTransient

noinst_PROGRAMS = $(EXAMPLE)

diagonalTransient_SOURCES = \
	 $(srcdir)/diagonalTransientMain.cpp \
	 $(top_srcdir)/../epetraext/example/model_evaluator/DiagonalTransient/EpetraExt_DiagonalTransientModel.cpp

diagonalTransient_DEPENDENCIES = $(common_dependencies)

diagonalTransient_LDADD = $(common_ldadd)


#
# Links
#

if ENABLE_RYTHMOS_STRATIMIKOS_AMESOS
LINK_DIAGONAL_AMESOS=diagonalTransient_amesos.exe
endif

if ENABLE_RYTHMOS_STRATIMIKOS_AZTECOO
LINK_DIAGONAL_AZTECOO=diagonalTransient_aztecoo.exe
if ENABLE_RYTHMOS_STRATIMIKOS_IFPACK
LINK_DIAGONAL_AZTECOO_IFPACK=diagonalTransient_aztecoo_ifpack.exe
endif
if ENABLE_RYTHMOS_STRATIMIKOS_ML
LINK_DIAGONAL_AZTECOO_ML=diagonalTransient_aztecoo_ml.exe
endif
endif

if ENABLE_RYTHMOS_STRATIMIKOS_BELOS
LINK_DIAGONAL_BELOS=diagonalTransient_belos.exe
if ENABLE_RYTHMOS_STRATIMIKOS_IFPACK
LINK_DIAGONAL_BELOS_IFPACK=diagonalTransient_belos_ifpack.exe
endif
if ENABLE_RYTHMOS_STRATIMIKOS_ML
LINK_DIAGONAL_BELOS_ML=diagonalTransient_belos_ml.exe
endif
endif

diagonalTransient_amesos.exe: diagonalTransient.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

diagonalTransient_aztecoo.exe: diagonalTransient.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

diagonalTransient_aztecoo_ifpack.exe: diagonalTransient.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

diagonalTransient_aztecoo_ml.exe: diagonalTransient.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

diagonalTransient_belos.exe: diagonalTransient.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

diagonalTransient_belos_ifpack.exe: diagonalTransient.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

diagonalTransient_belos_ml.exe: diagonalTransient.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

DIAGONAL_EXE_LINKS = \
  $(LINK_DIAGONAL_AMESOS) \
  $(LINK_DIAGONAL_AZTECOO) \
  $(LINK_DIAGONAL_AZTECOO_IFPACK) \
  $(LINK_DIAGONAL_AZTECOO_ML) \
  $(LINK_DIAGONAL_BELOS) \
  $(LINK_DIAGONAL_BELOS_IFPACK) \
  $(LINK_DIAGONAL_BELOS_ML)

#
# XML Files
#

XML_FILES = diagonalTransient_aztecoo.xml

diagonalTransient_aztecoo.xml: $(srcdir)/_diagonalTransient_aztecoo.xml
	cp $(srcdir)/_diagonalTransient_aztecoo.xml diagonalTransient_aztecoo.xml

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


EXTRA_DIST = _diagonalTransient_aztecoo.xml
