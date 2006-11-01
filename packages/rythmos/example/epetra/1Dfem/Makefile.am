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

#
# Common stuff
#

include $(top_builddir)/Makefile.export.rythmos
include $(top_builddir)/../nox/Makefile.export.nox

_EXPORT_INCLUDES = $(RYTHMOS_INCLUDES) $(STRATIMIKOS_INCLUDES) $(EPETRAEXT_INCLUDES) $(NOX_INCLUDES)
_EXPORT_LIBS = $(RYTHMOS_LIBS) $(STRATIMIKOS_LIBS) $(EPETRAEXT_LIBS) $(NOX_LIBS)


if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(_EXPORT_INCLUDES))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(_EXPORT_LIBS))
else
EXPORT_INCLUDES = $(_EXPORT_INCLUDES)
EXPORT_LIBS = $(_EXPORT_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

EXEEXT = .exe


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


