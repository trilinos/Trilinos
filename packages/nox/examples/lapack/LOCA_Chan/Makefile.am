## $Id$ 
## $Source$ 
#@HEADER
# ************************************************************************
# 
#                  LOCA Continuation Algorithm Package
#                 Copyright (2005) Sandia Corporation
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
# Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
# (rppawlo@sandia.gov), Sandia National Laboratories.
# 
# ************************************************************************
#@HEADER

# Include NOX export makefile
include $(top_builddir)/Makefile.export.nox

if USING_GNUMAKE
EXAMPLE_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(LOCA_LIBS))
EXAMPLE_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(LOCA_INCLUDES))
else
EXAMPLE_LIBS     = $(LOCA_LIBS)
EXAMPLE_INCLUDES = $(LOCA_INCLUDES)
endif

AM_CPPFLAGS = $(EXAMPLE_INCLUDES)

EXEEXT = .exe

CHANCONT = chanContinuation
CHANTPCONT = chanTPContinuation
CHANHOMOTOPY = chanHomotopy
CHANCONTCONSTRAINTS = chanContinuation_Constraints

if BUILD_LIBLOCAMF
CHANMULTICONT = chanMultiContinuation
endif

if BUILD_LOCA_LAPACK_EXAMPLES
noinst_PROGRAMS = $(CHANCONT) $(CHANTPCONT) $(CHANHOMOTOPY) \
	          $(CHANMULTICONT) $(CHANCONTCONSTRAINTS)
endif

chanContinuation_SOURCES = $(srcdir)/ChanContinuation.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanContinuation_DEPENDENCIES = \
        $(top_builddir)/src/libnox.a \
	$(top_builddir)/src-lapack/libnoxlapack.a \
	$(top_builddir)/src-loca/src/libloca.a \
        $(top_builddir)/src-loca/src-lapack/liblocalapack.a

chanContinuation_LDADD = $(EXAMPLE_LIBS)

chanTPContinuation_SOURCES = $(srcdir)/ChanTPContinuation.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanTPContinuation_DEPENDENCIES = $(chanContinuation_DEPENDENCIES)

chanTPContinuation_LDADD = $(EXAMPLE_LIBS)

chanHomotopy_SOURCES = $(srcdir)/ChanHomotopy.C \
		$(srcdir)/ChanProblemInterface.H \
		$(srcdir)/ChanProblemInterface.C

chanHomotopy_DEPENDENCIES = $(chanContinuation_DEPENDENCIES)

chanHomotopy_LDADD = $(EXAMPLE_LIBS)

chanMultiContinuation_SOURCES = \
	$(srcdir)/ChanMultiContinuation.C \
	$(srcdir)/ChanProblemInterface.H \
	$(srcdir)/ChanProblemInterface.C

chanMultiContinuation_DEPENDENCIES = $(chanContinuation_DEPENDENCIES) \
                                 $(top_builddir)/src-loca/src-mf/liblocamf.a 

chanMultiContinuation_LDADD = $(EXAMPLE_LIBS)

chanContinuation_Constraints_SOURCES = \
	$(srcdir)/ChanContinuation_Constraints.C \
	$(srcdir)/ChanConstraint.H \
	$(srcdir)/ChanConstraint.C \
	$(srcdir)/ChanProblemInterface.H \
	$(srcdir)/ChanProblemInterface.C

chanContinuation_Constraints_DEPENDENCIES = \
	$(chanContinuation_DEPENDENCIES)

chanContinuation_Constraints_LDADD = $(EXAMPLE_LIBS)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Files to be deleted by 'make clean'
# ------------------------------------------------------------------------

MOSTLYCLEANFILES = $(CXXREPOS)

