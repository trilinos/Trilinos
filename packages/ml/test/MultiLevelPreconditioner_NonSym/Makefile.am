#@HEADER
# ************************************************************************
# 
#                 ML
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

noinst_PROGRAMS = MultiLevelPreconditioner_NonSym 

if ML_USING_AMESOS
# Amesos require Teuchos
amesos_lib = $(top_builddir)/../amesos/src/libamesos.a
else
amesos_lib =
endif
if ML_USING_ANASAZI
# Anasazi require Teuchos
anasazi_lib = $(top_builddir)/../anasazi/src/libanasazi.a
else
anasazi_lib = 
endif
if ML_USING_TEUCHOS 
teuchos_lib = $(top_builddir)/../teuchos/src/libteuchos.a
else
teuchos_lib =
endif
if ML_USING_TRIUTILS
triutils_lib = $(top_builddir)/../triutils/src/libtriutils.a
else
triutils_lib =
endif
if ML_USING_IFPACK
ifpack_lib = $(top_builddir)/../ifpack/src/libifpack.a
else
ifpack_lib =
endif

MultiLevelPreconditioner_NonSym_SOURCES = $(srcdir)/MultiLevelPreconditioner_NonSym.cpp

MultiLevelPreconditioner_NonSym_DEPENDENCIES = \
	$(top_builddir)/src/libml.a \
	$(top_builddir)/../epetra/src/libepetra.a \
    ${triutils_lib} \
    ${amesos_lib} \
    ${anasazi_lib} \
    ${teuchos_lib} 

ml_common_include_locations = -I$(top_srcdir)/src/Include          \
                   -I$(top_srcdir)/src/Coarsen                  \
                   -I$(top_srcdir)/src/Comm                  \
                   -I$(top_srcdir)/src/FEGrid                  \
                   -I$(top_srcdir)/src/Include                  \
                   -I$(top_srcdir)/src/Krylov                  \
                   -I$(top_srcdir)/src/Main                  \
                   -I$(top_srcdir)/src/Operator                  \
                   -I$(top_srcdir)/src/Smoother                  \
                   -I$(top_srcdir)/src/Utils                  \
                   -I$(top_srcdir)/../aztecoo/src            \
                   -I$(top_builddir)/../aztecoo/src \
                   -I$(top_srcdir)/../epetra/src            \
                   -I$(top_builddir)/../epetra/src \
                   -I$(top_srcdir)/../teuchos/src            \
                   -I$(top_builddir)/../teuchos/src \
                   -I$(top_srcdir)/../triutils/src            \
                   -I$(top_builddir)/../triutils/src
#
ml_libs = $(top_builddir)/src/libml.a        \
                $(top_builddir)/../aztecoo/src/libaztecoo.a        \
                ${ifpack_lib} \
                ${amesos_lib} \
                ${anasazi_lib} \
                ${teuchos_lib} \
                ${triutils_lib} \
                $(top_builddir)/../epetra/src/libepetra.a \
                $(LIBS) \
                $(LAPACK_LIBS)                            \
                $(BLAS_LIBS)                            \
                $(FLIBS)

mllink = $(CXX) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

MultiLevelPreconditioner_NonSym_CPPFLAGS = ${ml_common_include_locations}
MultiLevelPreconditioner_NonSym_CXXFLAGS = ${ml_common_include_locations}
MultiLevelPreconditioner_NonSym_LDADD = ${ml_libs}
MultiLevelPreconditioner_NonSym_LINK = ${mllink}
