#@HEADER
# ************************************************************************
# 
#            Trilinos: An Object-Oriented Solver Framework
#                Copyright (2002) Sandia Corporation
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

noinst_PROGRAMS = AztecOO_msr adapt_msr AztecOO_vbr

AztecOO_msr_SOURCES = $(srcdir)/cxx_main.cpp

#AztecOO_msr_DEPENDENCIES = $(top_builddir)/src/libaztecoo.a \
#	$(top_builddir)/../epetra/src/libepetra.a\
#	$(top_builddir)/../triutils/src/libtriutils.a

AztecOO_msr_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/../epetra/src \
	-I$(top_srcdir)/../triutils/src -I$(top_builddir)/../epetra/src \
	-I$(top_builddir)/../triutils/src

AztecOO_msr_LDADD = $(top_builddir)/src/libaztecoo.a\
	$(top_builddir)/../triutils/src/libtriutils.a\
	$(top_builddir)/../epetra/src/libepetra.a\
	$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

AztecOO_vbr_SOURCES = $(srcdir)/vbr_main.cpp

#AztecOO_vbr_DEPENDENCIES = $(top_builddir)/src/libaztecoo.a \
#	$(top_builddir)/../epetra/src/libepetra.a\
#	$(top_builddir)/../triutils/src/libtriutils.a

AztecOO_vbr_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/../epetra/src \
	-I$(top_srcdir)/../triutils/src -I$(top_builddir)/../epetra/src \
	-I$(top_builddir)/../triutils/src

AztecOO_vbr_LDADD = $(top_builddir)/src/libaztecoo.a\
	$(top_builddir)/../triutils/src/libtriutils.a\
	$(top_builddir)/../epetra/src/libepetra.a\
	$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

adapt_msr_SOURCES = $(srcdir)/adapt_main.cpp

#adapt_msr_DEPENDENCIES = $(top_builddir)/src/libaztecoo.a \
#	$(top_builddir)/../epetra/src/libepetra.a\
#	$(top_builddir)/../triutils/src/libtriutils.a

adapt_msr_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/../epetra/src \
	-I$(top_srcdir)/../triutils/src -I$(top_builddir)/../epetra/src \
	-I$(top_builddir)/../triutils/src

adapt_msr_LDADD = $(top_builddir)/src/libaztecoo.a\
	$(top_builddir)/../triutils/src/libtriutils.a\
	$(top_builddir)/../epetra/src/libepetra.a\
	$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)


