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

all-local: make-links
	@echo ""
	@echo "rythmos/example/basicExample built successfully."
	@echo ""

#
# Actual programs
#

EXAMPLE = basicExample

noinst_PROGRAMS = $(EXAMPLE)

basicExample_SOURCES = \
	 $(srcdir)/cxx_main.cpp \
	 $(srcdir)/ExampleApplication.cpp

basicExample_DEPENDENCIES = $(common_dependencies)

basicExample_LDADD = $(common_ldadd)

EXTRA_DIST = $(srcdir)/ExampleApplication.hpp

#
# Links
#

if ENABLE_RYTHMOS_STRATIMIKOS_AMESOS
LINK_DIAGONAL_AMESOS=basicExample_amesos.exe
endif

if ENABLE_RYTHMOS_STRATIMIKOS_AZTECOO
LINK_DIAGONAL_AZTECOO=basicExample_aztecoo.exe
if ENABLE_RYTHMOS_STRATIMIKOS_IFPACK
LINK_DIAGONAL_AZTECOO_IFPACK=basicExample_aztecoo_ifpack.exe
endif
if ENABLE_RYTHMOS_STRATIMIKOS_ML
LINK_DIAGONAL_AZTECOO_ML=basicExample_aztecoo_ml.exe
endif
endif

if ENABLE_RYTHMOS_STRATIMIKOS_BELOS
LINK_DIAGONAL_BELOS=basicExample_belos.exe
if ENABLE_RYTHMOS_STRATIMIKOS_IFPACK
LINK_DIAGONAL_BELOS_IFPACK=basicExample_belos_ifpack.exe
endif
if ENABLE_RYTHMOS_STRATIMIKOS_ML
LINK_DIAGONAL_BELOS_ML=basicExample_belos_ml.exe
endif
endif

basicExample_amesos.exe:
	ln -s ./basicExample.exe basicExample_amesos.exe

basicExample_aztecoo.exe:
	ln -s ./basicExample.exe basicExample_aztecoo.exe

basicExample_aztecoo_ifpack.exe:
	ln -s ./basicExample.exe basicExample_aztecoo_ifpack.exe

basicExample_aztecoo_ml.exe:
	ln -s ./basicExample.exe basicExample_aztecoo_ml.exe

basicExample_belos.exe:
	ln -s ./basicExample.exe basicExample_belos.exe

basicExample_belos_ifpack.exe:
	ln -s ./basicExample.exe basicExample_belos_ifpack.exe

basicExample_belos_ml.exe:
	ln -s ./basicExample.exe basicExample_belos_ml.exe

make-links: \
  $(LINK_DIAGONAL_AMESOS) \
  $(LINK_DIAGONAL_AZTECOO) \
  $(LINK_DIAGONAL_AZTECOO_IFPACK) \
  $(LINK_DIAGONAL_AZTECOO_ML) \
  $(LINK_DIAGONAL_BELOS) \
  $(LINK_DIAGONAL_BELOS_IFPACK) \
  $(LINK_DIAGONAL_BELOS_ML)


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

# This is already added part of RYTHMOS_LIBS and therefore automake does not need to add this again!
LIBS =
