## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#            NOX: An Object-Oriented Nonlinear Solver Package
#                 Copyright (2002) Sandia Corporation
# 
#            LOCA: Library of Continuation Algorithms Package
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
# 
# Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
# Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
# ************************************************************************
#  CVS Information
#  $Source$
#  $Author$
#  $Date$
#  $Revision$
# ************************************************************************
# @HEADER

# Include NOX export makefile
include $(top_builddir)/Makefile.export.nox

if USING_GNUMAKE
EXAMPLE_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(NOX_LIBS))
EXAMPLE_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(NOX_INCLUDES))
else
EXAMPLE_LIBS     = $(NOX_LIBS)
EXAMPLE_INCLUDES = $(NOX_INCLUDES)
endif

AM_CPPFLAGS = $(EXAMPLE_INCLUDES)

EXEEXT = .exe

if BUILD_LAPACK_EXAMPLES
noinst_PROGRAMS = rosenbrock broyden $(PRERELEASE_EXAMPLES)
endif

# Rosenbrock Example
rosenbrock_SOURCES = $(srcdir)/Rosenbrock.C
rosenbrock_DEPENDENCIES = $(top_builddir)/src-lapack/libnoxlapack.a \
                          $(top_builddir)/src/libnox.a 
rosenbrock_LDADD = $(EXAMPLE_LIBS)

# Broyden Example
broyden_SOURCES = $(srcdir)/Broyden.C
broyden_DEPENDENCIES = $(rosenbrock_DEPENDENCIES)
broyden_LDADD = $(EXAMPLE_LIBS)

# Prerelease Linear Example
if BUILD_PRERELEASE
PRERELEASE_EXAMPLES = linear
endif
nodist_linear_SOURCES = $(srcdir)/Linear.C
linear_DEPENDENCIES = $(top_builddir)/src-lapack/libnoxlapack.a \
                      $(top_builddir)/src/libnox.a 
linear_LDADD = $(EXAMPLE_LIBS)

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libnoxexample_a_AR = $(ALTERNATE_AR)

else

libnoxexample_a_AR = $(AR) cru

endif

MAINTAINERCLEANFILES = Makefile.in 

MOSTLYCLEANFILES = $(CXXREPOS)


