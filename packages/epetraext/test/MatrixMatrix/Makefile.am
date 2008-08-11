#@HEADER
# ************************************************************************
# 
#     EpetraExt: Epetra Extended - Linear Algebra Services Package 
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

noinst_PROGRAMS = copyfiles MatrixMatrix_test

copyfiles_SOURCES=$(top_srcdir)/test/MatrixMatrix/src_infiles

MatrixMatrix_test_SOURCES = \
  $(srcdir)/cxx_main.cpp \
  $(top_srcdir)/test/epetra_test_err.h

MatrixMatrix_test_DEPENDENCIES = \
  $(top_builddir)/src/libepetraext.a \
  $(top_builddir)/../epetra/src/libepetra.a

include $(top_builddir)/Makefile.export.epetraext

if USING_GNUMAKE
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(EPETRAEXT_LIBS))
EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(EPETRAEXT_INCLUDES))
else
EXPORT_LIBS = $(EPETRAEXT_LIBS)
EXPORT_INC_PATH = $(EPETRAEXT_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INC_PATH)

MatrixMatrix_test_LDADD = $(EXPORT_LIBS)

EXTRA_DIST = \
	src_infileAB \
	src_infileATBT \
	src_infileAB2 \
	src_infiles \
	src_infileABT \
	src_infileAB3 \
	src_infileATB \
	src_infileATB2 \
	src_C.mtx \
	src_C4x4.mtx \
	src_C4x12x12x4.mtx \
	src_C4x6.mtx \
	src_C6x4.mtx \
	src_C6x6.mtx \
	src_mat6x4.mtx \
	src_mat6x6.mtx \
	src_mat12x4.mtx \
	src_mat4x12.mtx \
	src_mat4x4.mtx \
	src_mat4x6.mtx \
	src_Y.mtx \
	src_YTC.mtx \
	src_Y_transp.mtx \
	src_roman \
	src_roman2 \
	src_romancase \
	src_romancase2 \
	src_cdt \
	src_cdt_case \
	src_cdt_d.mtx \
	src_cdt_m.mtx \
	src_cdt_tce.mtx

copyfiles$(EXEEXT): $(copyfiles_SOURCES)
	cp $(top_srcdir)/test/MatrixMatrix/src_infileAB $(top_builddir)/test/MatrixMatrix/infileAB
	cp $(top_srcdir)/test/MatrixMatrix/src_infileATBT $(top_builddir)/test/MatrixMatrix/infileATBT
	cp $(top_srcdir)/test/MatrixMatrix/src_infileAB2 $(top_builddir)/test/MatrixMatrix/infileAB2
	cp $(top_srcdir)/test/MatrixMatrix/src_infiles $(top_builddir)/test/MatrixMatrix/infiles
	cp $(top_srcdir)/test/MatrixMatrix/src_infileABT $(top_builddir)/test/MatrixMatrix/infileABT
	cp $(top_srcdir)/test/MatrixMatrix/src_infileAB3 $(top_builddir)/test/MatrixMatrix/infileAB3
	cp $(top_srcdir)/test/MatrixMatrix/src_infileATB $(top_builddir)/test/MatrixMatrix/infileATB
	cp $(top_srcdir)/test/MatrixMatrix/src_infileATB2 $(top_builddir)/test/MatrixMatrix/infileATB2
	cp $(top_srcdir)/test/MatrixMatrix/src_C.mtx $(top_builddir)/test/MatrixMatrix/C.mtx
	cp $(top_srcdir)/test/MatrixMatrix/src_C4x4.mtx $(top_builddir)/test/MatrixMatrix/C4x4.mtx
	cp $(top_srcdir)/test/MatrixMatrix/src_C4x12x12x4.mtx $(top_builddir)/test/MatrixMatrix/C4x12x12x4.mtx
	cp $(top_srcdir)/test/MatrixMatrix/src_C4x6.mtx $(top_builddir)/test/MatrixMatrix/C4x6.mtx
	cp $(top_srcdir)/test/MatrixMatrix/src_C6x4.mtx $(top_builddir)/test/MatrixMatrix/C6x4.mtx
	cp $(top_srcdir)/test/MatrixMatrix/src_C6x6.mtx $(top_builddir)/test/MatrixMatrix/C6x6.mtx
	cp $(top_srcdir)/test/MatrixMatrix/src_mat6x4.mtx $(top_builddir)/test/MatrixMatrix/mat6x4.mtx
	cp $(top_srcdir)/test/MatrixMatrix/src_mat6x6.mtx $(top_builddir)/test/MatrixMatrix/mat6x6.mtx
	cp $(top_srcdir)/test/MatrixMatrix/src_mat12x4.mtx $(top_builddir)/test/MatrixMatrix/mat12x4.mtx
	cp $(top_srcdir)/test/MatrixMatrix/src_mat4x12.mtx $(top_builddir)/test/MatrixMatrix/mat4x12.mtx
	cp $(top_srcdir)/test/MatrixMatrix/src_mat4x4.mtx $(top_builddir)/test/MatrixMatrix/mat4x4.mtx
	cp $(top_srcdir)/test/MatrixMatrix/src_mat4x6.mtx $(top_builddir)/test/MatrixMatrix/mat4x6.mtx
	cp $(top_srcdir)/test/MatrixMatrix/src_Y.mtx $(top_builddir)/test/MatrixMatrix/Y.mtx
	cp $(top_srcdir)/test/MatrixMatrix/src_YTC.mtx $(top_builddir)/test/MatrixMatrix/YTC.mtx
	cp $(top_srcdir)/test/MatrixMatrix/src_Y_transp.mtx $(top_builddir)/test/MatrixMatrix/Y_transp.mtx
	cp $(top_srcdir)/test/MatrixMatrix/src_roman $(top_builddir)/test/MatrixMatrix/roman
	cp $(top_srcdir)/test/MatrixMatrix/src_roman2 $(top_builddir)/test/MatrixMatrix/roman2
	cp $(top_srcdir)/test/MatrixMatrix/src_romancase $(top_builddir)/test/MatrixMatrix/romancase
	cp $(top_srcdir)/test/MatrixMatrix/src_romancase2 $(top_builddir)/test/MatrixMatrix/romancase2
	cp $(top_srcdir)/test/MatrixMatrix/src_cdt $(top_builddir)/test/MatrixMatrix/cdt
	cp $(top_srcdir)/test/MatrixMatrix/src_cdt_case $(top_builddir)/test/MatrixMatrix/cdt_case
	cp $(top_srcdir)/test/MatrixMatrix/src_cdt_d.mtx $(top_builddir)/test/MatrixMatrix/cdt_d.mtx
	cp $(top_srcdir)/test/MatrixMatrix/src_cdt_m.mtx $(top_builddir)/test/MatrixMatrix/cdt_m.mtx
	cp $(top_srcdir)/test/MatrixMatrix/src_cdt_tce.mtx $(top_builddir)/test/MatrixMatrix/cdt_tce.mtx

