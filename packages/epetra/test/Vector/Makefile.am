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

noinst_PROGRAMS = Vector_test

Vector_test_SOURCES = $(srcdir)/cxx_main.cpp $(srcdir)/BuildTestProblems.cpp \
                  $(srcdir)/ExecuteTestProblems.cpp \
		  $(srcdir)/ExecuteTestProblems.h $(srcdir)/BuildTestProblems.h

Vector_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

Vector_test_CXXFLAGS = -I$(top_srcdir)/src

Vector_test_LDADD = $(top_builddir)/src/libepetra.a $(BLAS_LIBS) \
	$(LIBS) $(FLIBS)
