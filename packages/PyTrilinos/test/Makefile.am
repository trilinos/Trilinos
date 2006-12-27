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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
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
ACLOCAL	= aclocal -I config

EXTRA_DIST = definition setpath.py \
	     testAmesos.py testAnasazi.py testAztecOO.py testEpetraExt_Coloring.py \
	     testEpetra_BlockMap.py testEpetra_Comm.py testEpetra_CompObject.py \
	     testEpetra_CrsGraph.py testEpetra_CrsMatrix.py testEpetra_ImportExport.py \
	     testEpetra_IntSerialDense.py testEpetra_IntVector.py \
	     testEpetra_LinearProblem.py testEpetra_LocalMap.py testEpetra_Map.py \
	     testEpetra_MapColoring.py testEpetra_MultiVector.py testEpetra_Object.py \
	     testEpetra_SerialDense.py testEpetra_Vector.py testGaleri_CrsMatrices.py \
	     testGaleri_Maps.py testIFPACK.py testML_Preconditioner.py \
	     testNOX_StatusTest.py testTeuchos.py testTeuchos_ParameterList.py \
	     testThyra.py testTriUtils.py

MAINTAINERCLEANFILES = Makefile.in

## #######################################################################
## Define the suite of tests
## #######################################################################
if HAVE_PYTRILINOS_TEUCHOS
TEUCHOS_TESTS = testTeuchos testTeuchos_ParameterList
endif

if HAVE_PYTRILINOS_THYRA
THYRA_TESTS = testThyra
endif

if HAVE_PYTRILINOS_EPETRA
EPETRA_TESTS = testEpetra_BlockMap testEpetra_Comm testEpetra_CompObject \
	       testEpetra_CrsGraph testEpetra_CrsMatrix testEpetra_ImportExport \
	       testEpetra_IntSerialDense testEpetra_IntVector \
	       testEpetra_LinearProblem testEpetra_LocalMap testEpetra_Map \
	       testEpetra_MapColoring testEpetra_MultiVector testEpetra_Object \
	       testEpetra_SerialDense testEpetra_Vector
endif

if HAVE_PYTRILINOS_TRIUTILS
TRIUTILS_TESTS = testTriUtils
endif

if HAVE_PYTRILINOS_EPETRAEXT
EPETRAEXT_TESTS = testEpetraExt_Coloring
endif

if HAVE_PYTRILINOS_AZTECOO
AZTECOO_TESTS = testAztecOO
endif

if HAVE_PYTRILINOS_GALERI
GALERI_TESTS = testGaleri_CrsMatrices testGaleri_Maps
endif

if HAVE_PYTRILINOS_AMESOS
AMESOS_TESTS = testAmesos
endif

if HAVE_PYTRILINOS_IFPACK
IFPACK_TESTS = testIFPACK
endif

if HAVE_PYTRILINOS_ANASAZI
ANASAZI_TESTS = testAnasazi
endif

if HAVE_PYTRILINOS_ML
ML_TESTS = testML_Preconditioner
endif

if HAVE_PYTRILINOS_NOX
NOX_TESTS = testNOX_StatusTest
endif

TEST_SUITE = $(TEUCHOS_TESTS)   $(THYRA_TESTS)   $(EPETRA_TESTS) $(TRIUTILS_TESTS) \
	     $(EPETRAEXT_TESTS) $(AZTECOO_TESTS) $(GALERI_TESTS) $(AMESOS_TESTS)   \
	     $(IFPACK_TESTS)    $(ANASAZI_TESTS) $(ML_TESTS)     $(NOX_TESTS)

TEST_PYSCRIPTS = $(TEST_SUITE:=.py) 

SETPATH = setpath.py

# The following line helps the test harness recover from build errors.
all:
	@echo ""
	@echo "Trilinos package PyTrilinos subdirectory test built successfully."
	@echo ""

all-local: $(TEST_SUITE)
	@if [ `$(PYTHON) -c "import os; print os.path.abspath('$(srcdir)')"` != `pwd` ]; then \
	  echo "copying $(SETPATH)"; /bin/cp $(srcdir)/$(SETPATH) . ; fi

test%:
	@if [ `$(PYTHON) -c "import os; print os.path.abspath('$(srcdir)')"` != `pwd` ]; then \
	  echo "copying $@.py"; /bin/cp $(srcdir)/$@.py . ; fi

run-tests: all-local
	@for test_script in $(TEST_PYSCRIPTS); do if [ $$test_script != '.py' ]; then \
	  echo; echo "Running $$test_script"; echo; $$test_script -t; fi; done

CLEANFILES = $(TEST_PYSCRIPTS) $(SETPATH) *.pyc *.dat
