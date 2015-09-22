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

all-local: $(TIMEDISCR_EXE_LINKS) $(TIMEDISCR_INPUT_FILES)
	@echo ""
	@echo "rythmos/example/timeDiscretizedBackwardEuler built successfully."
	@echo ""


#
# Actual programs
#

EXAMPLE = timeDiscretizedBackwardEuler

noinst_PROGRAMS = $(EXAMPLE)

timeDiscretizedBackwardEuler_SOURCES = \
  $(srcdir)/timeDiscretizedBackwardEulerMain.cpp\
  $(top_srcdir)/../epetraext/example/model_evaluator/DiagonalTransient/EpetraExt_DiagonalTransientModel.cpp

timeDiscretizedBackwardEuler_DEPENDENCIES = $(common_dependencies)

timeDiscretizedBackwardEuler_LDADD = $(common_ldadd)


#
# Links
#

if ENABLE_RYTHMOS_STRATIMIKOS_AMESOS
LINK_TIMEDISCR_AMESOS=timeDiscretizedBackwardEuler_amesos.exe
endif

if ENABLE_RYTHMOS_STRATIMIKOS_AZTECOO
LINK_TIMEDISCR_AZTECOO=timeDiscretizedBackwardEuler_aztecoo.exe
if ENABLE_RYTHMOS_STRATIMIKOS_IFPACK
LINK_TIMEDISCR_AZTECOO_IFPACK=timeDiscretizedBackwardEuler_aztecoo_ifpack.exe
endif
if ENABLE_RYTHMOS_STRATIMIKOS_ML
LINK_TIMEDISCR_AZTECOO_ML=timeDiscretizedBackwardEuler_aztecoo_ml.exe
endif
endif

if ENABLE_RYTHMOS_STRATIMIKOS_BELOS
LINK_TIMEDISCR_BELOS=timeDiscretizedBackwardEuler_belos.exe
if ENABLE_RYTHMOS_STRATIMIKOS_IFPACK
LINK_TIMEDISCR_BELOS_IFPACK=timeDiscretizedBackwardEuler_belos_ifpack.exe
endif
if ENABLE_RYTHMOS_STRATIMIKOS_ML
LINK_TIMEDISCR_BELOS_ML=timeDiscretizedBackwardEuler_belos_ml.exe
endif
endif

timeDiscretizedBackwardEuler_amesos.exe: timeDiscretizedBackwardEuler.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

timeDiscretizedBackwardEuler_aztecoo.exe: timeDiscretizedBackwardEuler.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

timeDiscretizedBackwardEuler_aztecoo_ifpack.exe: timeDiscretizedBackwardEuler.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

timeDiscretizedBackwardEuler_aztecoo_ml.exe: timeDiscretizedBackwardEuler.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

timeDiscretizedBackwardEuler_belos.exe: timeDiscretizedBackwardEuler.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

timeDiscretizedBackwardEuler_belos_ifpack.exe: timeDiscretizedBackwardEuler.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

timeDiscretizedBackwardEuler_belos_ml.exe: timeDiscretizedBackwardEuler.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

TIMEDISCR_EXE_LINKS = \
  $(LINK_TIMEDISCR_AMESOS) \
  $(LINK_TIMEDISCR_AZTECOO) \
  $(LINK_TIMEDISCR_AZTECOO_IFPACK) \
  $(LINK_TIMEDISCR_AZTECOO_ML) \
  $(LINK_TIMEDISCR_BELOS) \
  $(LINK_TIMEDISCR_BELOS_IFPACK) \
  $(LINK_TIMEDISCR_BELOS_ML)


#
# XML input files
#

TIMEDISCR_INPUT_FILES = \
  timeDiscretizedBackwardEuler1.xml

timeDiscretizedBackwardEuler1.xml: $(srcdir)/_timeDiscretizedBackwardEuler1.xml
	cp $< $@


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


#
# All extra files that have to be added to the tarball!
#

EXTRA_DIST = \
  _timeDiscretizedBackwardEuler1.xml
