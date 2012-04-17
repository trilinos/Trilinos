#@HEADER
# ************************************************************************
# 
#               Epetra: Linear Algebra Services Package 
#                 Copyright (2001) Sandia Corporation
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

noinst_PROGRAMS = AztecOO

AztecOO_SOURCES = $(srcdir)/cxx_main.cpp

AztecOO_DEPENDENCIES = $(top_builddir)/src/libaztecoo.a

include $(top_builddir)/Makefile.export.aztecoo

if USING_GNUMAKE
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(AZTECOO_LIBS))
EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(AZTECOO_INCLUDES))
else
EXPORT_LIBS = $(AZTECOO_LIBS)
EXPORT_INC_PATH = $(AZTECOO_INCLUDES)
endif

##AztecOO_CXXFLAGS = -I$(top_srcdir)/src 

AM_CPPFLAGS = $(EXPORT_INC_PATH)

#Unsure why this test need to link to the BLAS on an IBM.
AztecOO_LDADD = $(EXPORT_LIBS)

