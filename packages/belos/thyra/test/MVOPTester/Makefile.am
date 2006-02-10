#@HEADER
# ************************************************************************
# 
#                 Belos: Block Linear Solvers Package
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

include $(top_builddir)/Makefile.export.belos

if USING_GNUMAKE
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(BELOS_TEST_LIBS))
EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(BELOS_TEST_INCLUDES))
else
EXPORT_LIBS = $(BELOS_TEST_LIBS)
EXPORT_INC_PATH = $(BELOS_TEST_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INC_PATH)

EXEEXT = .exe

#np# List all executable names here
noinst_PROGRAMS = MVOPTesterThyra_test

#np# One SOURCES primary is needed for EACH exe.  (Notice for all of the following 
#np# primaries that the pattern is "(Exe_name)_(PRIMARY_NAME) = ").  List all source 
#np# files needed to compile the exe.
MVOPTesterThyra_test_SOURCES = $(srcdir)/cxx_main.cpp

#np# One DEPENDENCIES primary is needed for each exe.  List any libraries that the 
#np# exe depends on.  (This primary is not absolutely essential.)
MVOPTesterThyra_test_DEPENDENCIES = \
	$(top_builddir)/src/libbelos.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a

#np# One CXXFLAGS primary can be defined for each exe.  Include any -I's necessary for 
#np# building the exe.
MVOPTesterThyra_test_CXXFLAGS = \
	-I$(top_srcdir)/thyra/src \
	-I$(top_srcdir)/src

#np# One LDADD primary can be defined for each exe.  Include any LDFLAGS necessary for 
#np# building the exe.
MVOPTesterThyra_test_LDADD = $(EXPORT_LIBS) 

# Dummy for new export makefiles
LIBS = 


