#@HEADER
# ************************************************************************
# 
#               Thyra: Interfaces and Support Code for the Interoperability of Abstract Numerical Algorithms
#                 Copyright (2004) Sandia Corporation
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
#@HEADER

EXEEXT = .exe

#
# Actual programs
#

SILLY_CG_SOLVE_SERIAL = sillyCgSolve_serial
SILLY_CG_SOLVE_MPI = sillyCgSolve_mpi
SILLY_POWER_METHOD_SERIAL = sillyPowerMethod_serial
EXAMPLE_IMPLICITLY_COMPOSED_LINEAR_OPERATORS = exampleImplicitlyComposedLinearOperators

noinst_PROGRAMS = $(SILLY_CG_SOLVE_SERIAL) $(SILLY_CG_SOLVE_MPI) $(SILLY_POWER_METHOD_SERIAL) \
  $(EXAMPLE_IMPLICITLY_COMPOSED_LINEAR_OPERATORS)

sillyCgSolve_serial_SOURCES = $(srcdir)/sillyCgSolve_serial.cpp
sillyCgSolve_serial_DEPENDENCIES = $(common_dependencies)
sillyCgSolve_serial_LDADD = $(common_ldadd)

sillyCgSolve_mpi_SOURCES = $(srcdir)/sillyCgSolve_mpi.cpp
sillyCgSolve_mpi_DEPENDENCIES = $(common_dependencies)
sillyCgSolve_mpi_LDADD = $(common_ldadd)

sillyPowerMethod_serial_SOURCES = $(srcdir)/sillyPowerMethod_serial.cpp
sillyPowerMethod_serial_DEPENDENCIES = $(common_dependencies)
sillyPowerMethod_serial_LDADD = $(common_ldadd)

exampleImplicitlyComposedLinearOperators_SOURCES = $(srcdir)/exampleImplicitlyComposedLinearOperators.cpp
exampleImplicitlyComposedLinearOperators_DEPENDENCIES = $(common_dependencies)
exampleImplicitlyComposedLinearOperators_LDADD = $(common_ldadd)


#
# Common stuff
#

common_dependencies = \
  $(top_builddir)/src/libthyra.a

include $(top_builddir)/Makefile.export.thyra

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(THYRA_INCLUDES))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(THYRA_LIBS))
else
EXPORT_INCLUDES = $(THYRA_INCLUDES)
EXPORT_LIBS = $(THYRA_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

common_ldadd = $(EXPORT_LIBS)

# We don't need to install these headers, but we need the files in the tarball.
EXTRA_DIST = \
  $(srcdir)/ExampleTridiagSerialLinearOp.hpp \
  $(srcdir)/ExampleTridiagSpmdLinearOp.hpp \
  $(srcdir)/sillierCgSolve.hpp \
  $(srcdir)/silliestCgSolve.hpp \
  $(srcdir)/sillyCgSolve.hpp \
  $(srcdir)/sillyModifiedGramSchmidt.hpp \
  $(srcdir)/sillyPowerMethod.hpp

# This is already added part of THYRA_LIBS and therefore automake does not need to add this again!
LIBS =
