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

include $(top_builddir)/Makefile.export.thyra

EXEEXT = .exe

#
# Actual programs
#

SILLY_CG_SOLVE_SERIAL = sillyCgSolve_serial
SILLY_CG_SOLVE_MPI = sillyCgSolve_mpi
SILLY_POWER_METHOD_SERIAL = sillyPowerMethod_serial

noinst_PROGRAMS = $(SILLY_CG_SOLVE_SERIAL) $(SILLY_CG_SOLVE_MPI) $(SILLY_POWER_METHOD_SERIAL)

sillyCgSolve_serial_SOURCES = $(srcdir)/sillyCgSolve_serial.cpp
sillyCgSolve_serial_DEPENDENCIES = $(common_dependencies)
sillyCgSolve_serial_LDADD = $(common_ldadd)

sillyCgSolve_mpi_SOURCES = $(srcdir)/sillyCgSolve_mpi.cpp
sillyCgSolve_mpi_DEPENDENCIES = $(common_dependencies)
sillyCgSolve_mpi_LDADD = $(common_ldadd)

sillyPowerMethod_serial_SOURCES = $(srcdir)/sillyPowerMethod_serial.cpp
sillyPowerMethod_serial_DEPENDENCIES = $(common_dependencies)
sillyPowerMethod_serial_LDADD = $(common_ldadd)

#
# Common stuff
#

common_dependencies = \
  $(top_builddir)/src/libthyra.a

AM_CPPFLAGS = $(shell perl $(top_srcdir)/config/strip_dup_incl_paths.pl $(THYRA_INCLUDES))

common_ldadd = $(shell perl $(top_srcdir)/config/strip_dup_libs.pl $(THYRA_LIBS))

# This is already added part of THYRA_LIBS and therefore automake does not need to add this again!
LIBS =
