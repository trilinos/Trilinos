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
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or
# Ray Tuminaro (rstumin@sandia.gov) 
# 
# ************************************************************************
#@HEADER

EXEEXT = .exe

# One of the input options for this example is the coarse direct solver KLU.
# KLU is part of amesos, which in turn needs teuchos.
if ML_USING_AMESOS
if ML_USING_TEUCHOS

EXTRA_DIST= 2d_Poisson.csh baseline1 baseline2 baseline3 baseline5 baseline6 baseline7 ml_inputfile1 ml_inputfile2 ml_inputfile3 ml_inputfile5 ml_inputfile6 ml_inputfile7

noinst_PROGRAMS = 2d_Poisson

2d_Poisson_SOURCES = $(top_srcdir)/examples/ml_example2d.c

2d_Poisson$(EXT): $(2d_Poisson_SOURCES)

#copy these files over to the build area
GetNecessaryFiles: $(EXTRA_DIST)
	cp -f $(srcdir)/2d_Poisson.csh .
	cp -f $(srcdir)/baseline1 .
	cp -f $(srcdir)/baseline2 .
	cp -f $(srcdir)/baseline3 .
	cp -f $(srcdir)/baseline5 .
	cp -f $(srcdir)/baseline6 .
	cp -f $(srcdir)/baseline7 .
	cp -f $(srcdir)/ml_inputfile1 .
	cp -f $(srcdir)/ml_inputfile2 .
	cp -f $(srcdir)/ml_inputfile3 .
	cp -f $(srcdir)/ml_inputfile5 .
	cp -f $(srcdir)/ml_inputfile6 .
	cp -f $(srcdir)/ml_inputfile7 .

#2d_Poisson$(REALEXT): $(2d_Poisson_SOURCES)

if ML_USING_TEUCHOS 
teuchos_lib = $(top_builddir)/../teuchos/src/libteuchos.a
else
teuchos_lib =
endif

if ML_USING_AMESOS
# Amesos require Teuchos
amesos_lib = $(top_builddir)/../amesos/src/libamesos.a \
             $(top_builddir)/../teuchos/src/libteuchos.a
else
amesos_lib =
endif

if ML_USING_IFPACK
ifpack_lib = $(top_builddir)/../ifpack/src/libifpack.a
else
ifpack_lib =
endif

2d_Poisson_DEPENDENCIES = \
	$(top_builddir)/src/libml.a \
	$(top_builddir)/../aztecoo/src/libaztecoo.a \
    ${amesos_lib} \
	$(top_builddir)/../epetra/src/libepetra.a \
    ${teuchos_lib} \
	GetNecessaryFiles

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
                   -I$(top_builddir)/../epetra/src          \
                   -I$(top_builddir)/../teuchos/src        \
                   -I$(top_builddir)/../amesos/src

ml_libs = $(top_builddir)/src/libml.a        \
                $(top_builddir)/../aztecoo/src/libaztecoo.a        \
                ${ifpack_lib} \
	            ${amesos_lib} \
                $(top_builddir)/../epetra/src/libepetra.a \
	            ${teuchos_lib} \
                $(LIBS) \
                $(LAPACK_LIBS)                            \
                $(BLAS_LIBS)                            \
                $(FLIBS)

mllink = $(CXX) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

2d_Poisson_CPPFLAGS = ${ml_common_include_locations} -DML_BENCHMARK
2d_Poisson_CXXFLAGS = ${ml_common_include_locations} -DML_BENCHMARK
2d_Poisson_LDADD = ${ml_libs}
2d_Poisson_LINK = ${mllink}

endif #if ML_USING_TEUCHOS
endif #if ML_USING_AMESOS
