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

if HAVE_EPETRAEXT_HDF5
if HAVE_EPETRAEXT_TEUCHOS
HDF5Files = HDF5_IO HDF5Converter HDF5_MATLAB
include $(top_builddir)/../teuchos/Makefile.export.teuchos
else
HDF5Files = 
TEUCHOS_LIBS =
TEUCHOS_INCLUDES =
endif
else
HDF5Files = 
TEUCHOS_LIBS =
TEUCHOS_INCLUDES =
endif

noinst_PROGRAMS = MatrixMarket_IO XML_IO $(HDF5Files) 

MatrixMarket_IO_SOURCES = $(srcdir)/MatrixMarket_IO.cpp
MatrixMarket_IO_DEPENDENCIES = $(top_builddir)/src/libepetraext.a \
				  $(top_builddir)/../epetra/src/libepetra.a
MatrixMarket_IO_CXXFLAGS = -I$(top_srcdir)/src/inout \
			      -I$(top_srcdir)/src \
			      -I$(top_builddir)/src \
			      -I$(top_srcdir)/../epetra/src \
			      -I$(top_builddir)/../epetra/src \
			      -I$(top_srcdir)/../triutils/src \
			      -I$(top_builddir)/../triutils/src
MatrixMarket_IO_LDADD = $(top_builddir)/src/libepetraext.a \
			   $(top_builddir)/../epetra/src/libepetra.a \
			   $(top_builddir)/../triutils/src/libtriutils.a \
			   $(LIBS) \
			   $(BLAS_LIBS) \
			   $(FLIBS)

XML_IO_SOURCES = $(srcdir)/XML_IO.cpp
XML_IO_DEPENDENCIES = $(top_builddir)/src/libepetraext.a \
		      $(top_builddir)/../epetra/src/libepetra.a
XML_IO_CXXFLAGS = -I$(top_srcdir)/src/inout \
		  -I$(top_srcdir)/src \
		  -I$(top_builddir)/src \
		  -I$(top_srcdir)/../epetra/src \
		  -I$(top_builddir)/../epetra/src \
		  -I$(top_srcdir)/../teuchos/src \
		  -I$(top_builddir)/../teuchos/src
XML_IO_LDADD = $(top_builddir)/src/libepetraext.a \
	       $(top_builddir)/../epetra/src/libepetra.a \
	       $(top_builddir)/../teuchos/src/libteuchos.a \
	       $(LIBS) \
	       $(BLAS_LIBS) \
	       $(FLIBS)

HDF5_IO_SOURCES = $(srcdir)/HDF5_IO.cpp
HDF5_IO_DEPENDENCIES = $(top_builddir)/src/libepetraext.a \
			  $(top_builddir)/../epetra/src/libepetra.a
HDF5_IO_CXXFLAGS = -I$(top_srcdir)/src/inout \
		      -I$(top_srcdir)/src \
		      -I$(top_builddir)/src \
		      -I$(top_srcdir)/../epetra/src \
		      -I$(top_builddir)/../epetra/src \
		      -I$(top_srcdir)/../triutils/src \
		      -I$(top_builddir)/../triutils/src \
		      $(TEUCHOS_INCLUDES) 
HDF5_IO_LDADD = $(top_builddir)/src/libepetraext.a \
		   $(top_builddir)/../epetra/src/libepetra.a \
		   $(top_builddir)/../triutils/src/libtriutils.a \
		   $(TEUCHOS_LIBS) \
		   $(LIBS) \
		   $(BLAS_LIBS) \
		   $(FLIBS)

HDF5_MATLAB_SOURCES = $(srcdir)/HDF5_MATLAB.cpp
HDF5_MATLAB_DEPENDENCIES = $(top_builddir)/src/libepetraext.a \
			   $(top_builddir)/../epetra/src/libepetra.a
HDF5_MATLAB_CXXFLAGS = -I$(top_srcdir)/src/inout \
		       -I$(top_srcdir)/src \
		       -I$(top_builddir)/src \
		       -I$(top_srcdir)/../epetra/src \
		       -I$(top_builddir)/../epetra/src \
		       -I$(top_srcdir)/../triutils/src \
		       -I$(top_builddir)/../triutils/src \
		       $(TEUCHOS_INCLUDES) 
HDF5_MATLAB_LDADD = $(top_builddir)/src/libepetraext.a \
		    $(top_builddir)/../epetra/src/libepetra.a \
		    $(top_builddir)/../triutils/src/libtriutils.a \
		    $(TEUCHOS_LIBS) \
		    $(LIBS) \
		    $(BLAS_LIBS) \
		    $(FLIBS)

HDF5Converter_SOURCES = $(srcdir)/HDF5Converter.cpp
HDF5Converter_DEPENDENCIES = $(top_builddir)/src/libepetraext.a \
			     $(top_builddir)/../epetra/src/libepetra.a
HDF5Converter_CXXFLAGS = -I$(top_srcdir)/src/inout \
			 -I$(top_srcdir)/src \
			 -I$(top_builddir)/src \
			 -I$(top_srcdir)/../epetra/src \
			 -I$(top_builddir)/../epetra/src \
			 -I$(top_srcdir)/../triutils/src \
			 -I$(top_builddir)/../triutils/src \
			 $(TEUCHOS_INCLUDES) 
HDF5Converter_LDADD = $(top_builddir)/src/libepetraext.a \
		      $(top_builddir)/../epetra/src/libepetra.a \
		      $(top_builddir)/../triutils/src/libtriutils.a \
		      $(TEUCHOS_LIBS) \
		      $(LIBS) \
		      $(BLAS_LIBS) \
		      $(FLIBS)
