#@HEADER
# ************************************************************************
# 
#                     New_Package Example Package
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

#np# This file is an example of a typical Makefile.am to build an example exe. 

EXEEXT = .exe

#np# List the names of all exe's in this directory here
noinst_PROGRAMS = example

#np# One SOURCES primary is needed for EACH exe listed above.  (Notice for all 
#np# of the  following primaries that the pattern is 
#np# "(Exe_name)_(PRIMARY_NAME) = ").
#np# List all source files needed to compile the exe.
example_SOURCES = $(srcdir)/cxx_main.cpp

#np# One DEPENDENCIES primary is needed for each exe listed above.  List any
#np# libraries that the exe depends on.  (This primary is not absolutely
#np# essential.)
example_DEPENDENCIES = $(top_builddir)/src/libepetra.a

#np# One CXXFLAGS primary can be defined for each exe listed above.  Include 
#np# any -I's necessary for building the exe.
example_CXXFLAGS = -I$(top_srcdir)/src

#np# One LDADD primary can be defined for each exe listed above.  Include any
#np# LDFLAGS  necessary for building the exe.
example_LDADD = $(top_builddir)/src/libepetra.a \
	$(LIBS) $(BLAS_LIBS) $(FLIBS)
