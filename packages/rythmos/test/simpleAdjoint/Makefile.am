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

all-local: $(SIMPLE_ADJOINT_XML_FILES) $(SIMPLE_ADJOINT_EXE_LINKS)
	@echo ""
	@echo "rythmos/test/simpleAdjoint built successfully."
	@echo ""


#
# Actual programs
#

EXAMPLE = simpleAdjoint

noinst_PROGRAMS = $(EXAMPLE)

simpleAdjoint_SOURCES = \
	 $(srcdir)/simpleAdjointMain.cpp \
	 $(top_srcdir)/../epetraext/example/model_evaluator/DiagonalTransient/EpetraExt_DiagonalTransientModel.cpp

simpleAdjoint_DEPENDENCIES = $(common_dependencies)

simpleAdjoint_LDADD = $(common_ldadd)


#
# XML Files
#

SIMPLE_ADJOINT_XML_FILES = simpleAdjoint_n4.xml

simpleAdjoint_n4.xml: _simpleAdjoint_n4.xml
	cp $(srcdir)/_simpleAdjoint_n4.xml simpleAdjoint_n4.xml


#
# Links
#

if ENABLE_RYTHMOS_STRATIMIKOS_AMESOS
LINK_SIMPLE_ADJOINT_AMESOS=simpleAdjoint_amesos.exe
endif

simpleAdjoint_amesos.exe: simpleAdjoint.exe
	if [ ! -x $@ ]; then ln -s $< $@; fi

SIMPLE_ADJOINT_EXE_LINKS = \
  $(LINK_SIMPLE_ADJOINT_AMESOS)

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


# Add files to tarball that are not explictly know by autotools
EXTRA_DIST = _simpleAdjoint_n4.xml
