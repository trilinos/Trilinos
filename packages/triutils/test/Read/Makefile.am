#@HEADER
# ************************************************************************
# 
#                 TriUtils: Trilinos Utilities Package
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

#
#

noinst_PROGRAMS = read_test Dummy


Dummy$(EXEEXT): $(srcdir)/Dummy
	cp  $(srcdir)/In_bcsstk01.rsa  bcsstk01.rsa
	cp  $(srcdir)/In_bcsstk01.mtx  bcsstk01.mtx
	cp  $(srcdir)/In_bcsstk01.triS  bcsstk01.triS

	cp  $(srcdir)/In_Diagonal.rua  Diagonal.rua
	cp  $(srcdir)/In_Diagonal.mtx  Diagonal.mtx
	cp  $(srcdir)/In_Diagonal.TimD  Diagonal.TimD

	cp  $(srcdir)/In_fs_183_4.rua  fs_183_4.rua
	cp  $(srcdir)/In_fs_183_4.mtx  fs_183_4.mtx
	cp  $(srcdir)/In_fs_183_4.triU  fs_183_4.triU

	cp  $(srcdir)/In_impcol_a.rua  impcol_a.rua
	cp  $(srcdir)/In_impcol_a.mtx  impcol_a.mtx
	cp  $(srcdir)/In_impcol_a.triU  impcol_a.triU




SRC_CXX =  $(srcdir)/TestTriutilReads.cpp

read_test_SOURCES = 	 $(SRC_CXX)

Dummy_SOURCES = Dummy

include $(top_builddir)/Makefile.export.triutils

if USING_GNUMAKE
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(TRIUTILS_LIBS))
EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(TRIUTILS_INCLUDES))
else
EXPORT_LIBS = $(TRIUTILS_LIBS)
EXPORT_INC_PATH = $(TRIUTILS_INCLUDES)
endif

read_test_LDADD = $(EXPORT_LIBS)
AM_CPPFLAGS = $(EXPORT_INC_PATH)

read_test_DEPENDENCIES = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/src/libtriutils.a 

read_test_CXXFLAGS = -I$(top_srcdir)/src

EXTRA_DIST = \
	$(srcdir)/In_bcsstk01.rsa $(srcdir)/In_bcsstk01.mtx \
	$(srcdir)/In_bcsstk01.triS $(srcdir)/In_fs_183_4.rua \
	$(srcdir)/In_fs_183_4.mtx $(srcdir)/In_fs_183_4.triU \
	$(srcdir)/In_Diagonal.rua $(srcdir)/In_Diagonal.mtx \
	$(srcdir)/In_Diagonal.TimD \
	$(srcdir)/In_impcol_a.rua $(srcdir)/In_impcol_a.mtx \
	$(srcdir)/In_impcol_a.triU 
