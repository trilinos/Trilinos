#@HEADER
# ************************************************************************
# 
#            Trilinos: An Object-Oriented Solver Framework
#                Copyright (2001) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#   
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#   
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# 
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
# 
# ************************************************************************
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = inout_example

inout_example_SOURCES = \
  $(srcdir)/cxx_main.cpp

inout_example_DEPENDENCIES = \
  $(top_builddir)/src/transform/libepetraext.a \
  $(top_builddir)/../epetra/src/libepetra.a

inout_example_CXXFLAGS = -I$(top_srcdir)/src/inout \
                         -I$(top_builddir)/src \
                         -I$(top_srcdir)/../epetra/src \
                         -I$(top_builddir)/../epetra/src \
                         -I$(top_srcdir)/../triutils/src \
                         -I$(top_builddir)/../triutils/src

inout_example_LDADD = \
  $(top_builddir)/src/transform/libepetraext.a \
  $(top_builddir)/../epetra/src/libepetra.a \
  $(top_builddir)/../epetra/src/libtriutils.a \
  $(BLAS_LIBS) \
  $(LIBS) \
  $(FLIBS)

