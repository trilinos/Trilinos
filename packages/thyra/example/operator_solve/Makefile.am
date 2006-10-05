#@HEADER
# ************************************************************************
# 
#               Thyra: Interfaces and Support Code for the Interoperability of Abstract Numerical Algorithms
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

#
# Actual programs
#

if HAVE_TEUCHOS_COMPLEX
SILLY_1D_FFT_SERIAL = silly1DFFT_serial
else
SILLY_1D_FFT_SERIAL =
endif

noinst_PROGRAMS = $(SILLY_1D_FFT_SERIAL)

silly1DFFT_serial_SOURCES = $(srcdir)/silly1DFFT_serial.cpp
silly1DFFT_serial_DEPENDENCIES = $(common_dependencies)
silly1DFFT_serial_LDADD = $(common_ldadd)

#
# Common stuff
#

common_dependencies = \
  $(top_builddir)/src/libthyra.a

include $(top_builddir)/Makefile.export.thyra

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(THYRA_INCLUDES))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(THYRA_LIBS))
else
EXPORT_INCLUDES = $(THYRA_INCLUDES)
EXPORT_LIBS = $(THYRA_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

common_ldadd = $(EXPORT_LIBS)

# We don't need to install these headers, but we need the files in the tarball.
EXTRA_DIST = \
	     $(srcdir)/serial_1D_FFT.hpp \
	     $(srcdir)/RealComplexFFTLinearOp.hpp \
	     $(srcdir)/ComplexFFTLinearOp.hpp \
	     $(srcdir)/Thyra_LinearOpWithSolveFactoryExamples.hpp

# This is already added part of THYRA_LIBS and therefore automake does not need to add this again!
LIBS =
