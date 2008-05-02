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

include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
TEST_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

if ML_USING_TEUCHOS
if ML_USING_AZTECOO
if ML_USING_EPETRA
if ML_USING_GALERI
if ML_USING_AMESOS

noinst_PROGRAMS = ml_Validate GetNecessaryFiles

GetNecessaryFiles_SOURCES= $(top_srcdir)/test/ValidateParams

#copy these files over to the build area
GetNecessaryFiles$(EXEEXT): $(GetNecessaryFiles_SOURCES)
	@if test ! -e ./input1.xml; then cp -f $(srcdir)/input1.xml .; fi
	@if test ! -e ./baseline1.txt; then cp -f $(srcdir)/baseline1.txt .; fi
	@if test ! -e ./evaluateCriteria; then cp -f $(srcdir)/evaluateCriteria .; fi

ml_Validate_SOURCES      = $(srcdir)/cxx_main.cpp
ml_Validate_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_Validate_CXXFLAGS     = ${TEST_INCLUDES}
ml_Validate_LDADD        = ${TEST_LIBS}

endif
endif
endif
endif
endif
