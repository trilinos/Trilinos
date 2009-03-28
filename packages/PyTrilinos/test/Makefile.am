# @HEADER
# ************************************************************************
# 
#                PyTrilinos: Python Interface to Trilinos
#                   Copyright (2005) Sandia Corporation
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
# Questions? Contact Bill Spotz (wfspotz@sandia.gov) 
# 
# ************************************************************************
# @HEADER

## #######################################################################
## Options to automake (rarely used - don't worry about it)
## #######################################################################
AUTOMAKE_OPTIONS = foreign

## #######################################################################
## Aclocal command (rarely used - don't worry about it)
## #######################################################################
ACLOCAL_AMFLAGS = -I config

EXTRA_DIST = definition setpath.src \
	     testAmesos.src testAnasazi.src testAztecOO.src \
	     testEpetraExt_Coloring.src testEpetra_BlockMap.src \
	     testEpetra_Comm.src testEpetra_CompObject.src \
	     testEpetra_CrsGraph.src testEpetra_CrsMatrix.src \
	     testEpetra_FEVector.src testEpetra_ImportExport.src \
	     testEpetra_IntSerialDense.src testEpetra_IntVector.src \
	     testEpetra_LinearProblem.src testEpetra_LocalMap.src \
	     testEpetra_Map.src testEpetra_MapColoring.src \
	     testEpetra_MultiVector.src testEpetra_Object.src \
	     testEpetra_SerialDense.src testEpetra_Util.src \
	     testEpetra_Vector.src testGaleri_CrsMatrices.src \
	     testGaleri_Maps.src testIFPACK.src testKomplex_LinearProblem.src \
	     testML_Space.src testML_CompObject.src testML_TimeObject.src \
	     testML_MultiVector.src testML_PyMatrix.src \
	     testML_Preconditioner.src testNOX_StatusTest.src testTeuchos.src \
	     testTeuchos_Comm.src testTeuchos_ParameterList.src \
	     testTeuchos_ScalarTraits.src testThyra.src testTriUtils.src

MAINTAINERCLEANFILES = Makefile.in

## #######################################################################
## Define the suite of tests
## #######################################################################
if HAVE_PYTRILINOS_TEUCHOS
TEUCHOS_TESTS = testTeuchos.py testTeuchos_Comm.py \
	        testTeuchos_ParameterList.py testTeuchos_ScalarTraits.py
endif

# if HAVE_PYTRILINOS_THYRA
# THYRA_TESTS = testThyra
# endif

if HAVE_PYTRILINOS_EPETRA
EPETRA_TESTS = testEpetra_BlockMap.py testEpetra_Comm.py \
	       testEpetra_CompObject.py testEpetra_CrsGraph.py \
	       testEpetra_CrsMatrix.py testEpetra_FEVector.py \
	       testEpetra_ImportExport.py testEpetra_IntSerialDense.py \
	       testEpetra_IntVector.py testEpetra_LinearProblem.py \
	       testEpetra_LocalMap.py testEpetra_Map.py \
	       testEpetra_MapColoring.py testEpetra_MultiVector.py \
	       testEpetra_Object.py testEpetra_SerialDense.py \
	       testEpetra_Util.py testEpetra_Vector.py
endif

if HAVE_PYTRILINOS_TRIUTILS
TRIUTILS_TESTS = testTriUtils.py
endif

if HAVE_PYTRILINOS_EPETRAEXT
EPETRAEXT_TESTS = testEpetraExt_Coloring.py
endif

if HAVE_PYTRILINOS_AZTECOO
AZTECOO_TESTS = testAztecOO.py
endif

if HAVE_PYTRILINOS_GALERI
GALERI_TESTS = testGaleri_CrsMatrices.py testGaleri_Maps.py
endif

if HAVE_PYTRILINOS_AMESOS
AMESOS_TESTS = testAmesos.py
endif

if HAVE_PYTRILINOS_IFPACK
IFPACK_TESTS = testIFPACK.py
endif

if HAVE_PYTRILINOS_KOMPLEX
KOMPLEX_TESTS = testKomplex_LinearProblem.py
endif

if HAVE_PYTRILINOS_ANASAZI
ANASAZI_TESTS = testAnasazi.py
endif

if HAVE_PYTRILINOS_ML
ML_TESTS = testML_Space.py testML_CompObject.py testML_TimeObject.py \
	   testML_MultiVector.py testML_PyMatrix.py testML_Preconditioner.py
endif

if HAVE_PYTRILINOS_NOX
NOX_TESTS = testNOX_StatusTest.py
endif

TEST_SUITE = $(TEUCHOS_TESTS)   $(THYRA_TESTS)   $(EPETRA_TESTS)  $(TRIUTILS_TESTS) \
	     $(EPETRAEXT_TESTS) $(AZTECOO_TESTS) $(GALERI_TESTS)  $(AMESOS_TESTS) \
	     $(IFPACK_TESTS)    $(KOMPLEX_TESTS) $(ANASAZI_TESTS) $(ML_TESTS) \
	     $(NOX_TESTS)

TESTUTIL   = testutil.py

# The following line helps the test harness recover from build errors.
all:

all-local: $(TESTUTIL) $(TEST_SUITE)

%.py: %.src
	@echo "copying $@"; /bin/cp $< $@

run-tests: all-local
	@for test_script in $(TEST_SUITE); do echo; echo "Running $$test_script"; \
	  echo; $$test_script -t; done

CLEANFILES = $(TEST_SUITE) $(TESTUTIL) *.pyc *.dat *.xml
