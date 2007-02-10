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
ACLOCAL	= aclocal -I config

EXTRA_DIST = setpath.py gre__115.rua \
	     exAmesos_Factory.py exAmesos_Solvers.py exAztecOO.py exAztecOO_Operator.py \
	     exAztecOO_RowMatrix.py exEpetra.py exEpetra_Comm.py \
	     exEpetra_CrsMatrix_Easy.py exEpetra_CrsMatrix_Efficient.py \
	     exEpetra_FECrsMatrix_Easy.py exEpetra_ImportExport.py \
	     exEpetraExt_IO_MatrixMarket.py exEpetraExt_IO_XML.py \
	     exEpetraExt_MatrixMatrix.py exGaleri.py exGaleri_ReadHB.py exIFPACK.py \
	     exMLAPI.py exMLAPI_Simple.py exMLAPI_AztecOO.py exMLAPI_Iterate.py \
	     exMLAPI_PyMatrix.py exMLAPI_Smoother.py exNOX_1Dfem.py exThyra.py

MAINTAINERCLEANFILES = Makefile.in

## #######################################################################
## Define the suite of examples
## #######################################################################

# Teuchos examples
if HAVE_PYTRILINOS_TEUCHOS
TEUCHOS_EXAMPLES = 
endif

# # Thyra examples
# if HAVE_PYTRILINOS_THYRA
# if HAVE_PYTRILINOS_EPETRA
# THYRA_EXAMPLES = exThyra
# endif
# endif

# Epetra examples
if HAVE_PYTRILINOS_EPETRA
EPETRA_EXAMPLES = exEpetra exEpetra_Comm exEpetra_CrsMatrix_Easy \
	          exEpetra_CrsMatrix_Efficient exEpetra_FECrsMatrix_Easy \
	          exEpetra_ImportExport
endif

# TriUtils examples
if HAVE_PYTRILINOS_TRIUTILS
TRIUTILS_EXAMPLES = 
endif

# EpetraExt examples
if HAVE_PYTRILINOS_EPETRAEXT
if HAVE_PYTRILINOS_EPETRA
EPETRAEXT_EXAMPLES = exEpetraExt_IO_MatrixMarket exEpetraExt_IO_XML \
	             exEpetraExt_MatrixMatrix
endif
endif

# AztecOO examples
if HAVE_PYTRILINOS_AZTECOO
if HAVE_PYTRILINOS_EPETRA
if HAVE_PYTRILINOS_GALERI
AZTECOO_EX1 = exAztecOO
endif
AZTECOO_EX2 = exAztecOO_Operator exAztecOO_RowMatrix
AZTECOO_EXAMPLES = $(AZTECOO_EX1) $(AZTECOO_EX2)
endif
endif

# Galeri examples
if HAVE_PYTRILINOS_GALERI
if HAVE_PYTRILINOS_EPETRA
GALERI_EXAMPLES = exGaleri exGaleri_ReadHB
endif
endif

# Amesos examples
if HAVE_PYTRILINOS_AMESOS
if HAVE_PYTRILINOS_EPETRA
AMESOS_EXAMPLES = exAmesos_Factory exAmesos_Solvers
endif
endif

# IFPACK examples
if HAVE_PYTRILINOS_IFPACK
if HAVE_PYTRILINOS_AZTECOO
if HAVE_PYTRILINOS_GALERI
if HAVE_PYTRILINOS_TRIUTILS
if HAVE_PYTRILINOS_EPETRA
IFPACK_EXAMPLES = exIFPACK
endif
endif
endif
endif
endif

# Anasazi examples
if HAVE_PYTRILINOS_ANASAZI
ANASAZI_EXAMPLES = 
endif

# ML examples
if HAVE_PYTRILINOS_ML
if HAVE_PYTRILINOS_EPETRA
ML_EX1 = exMLAPI exMLAPI_Iterate exMLAPI_PyMatrix exMLAPI_Smoother
if HAVE_PYTRILINOS_AZTECOO
ML_EX2 = exMLAPI_AztecOO
if HAVE_PYTRILINOS_GALERI
ML_EX3 = exMLAPI_Simple exMLAPI_AztecOO
endif
endif
ML_EXAMPLES = $(ML_EX1) $(ML_EX2) $(ML_EX3)
endif
endif

# NOX examples
if HAVE_PYTRILINOS_NOX
if HAVE_PYTRILINOS_NOX_EPETRA
NOX_EXAMPLES = exNOX_1Dfem
endif
endif

# Collect the package examples into a suite of examples
EXAMPLE_SUITE = $(TEUCHOS_EXAMPLES)  $(THYRA_EXAMPLES)     $(EPETRA_EXAMPLES) \
	        $(TRIUTILS_EXAMPLES) $(EPETRAEXT_EXAMPLES) $(AZTECOO_EXAMPLES) \
	        $(GALERI_EXAMPLES)   $(AMESOS_EXAMPLES)    $(IFPACK_EXAMPLES) \
	        $(ANASAZI_EXAMPLES)  $(ML_EXAMPLES)        $(NOX_EXAMPLES)

EXAMPLE_PYSCRIPTS = $(EXAMPLE_SUITE:=.py) 

SHOULD_COPY = should_copy.py
SETPATH     = setpath.py
GRE__115    = gre__115.rua

# The following line helps the test harness recover from build errors.
all:
	@echo ""
	@echo "Trilinos package PyTrilinos subdirectory example built successfully."
	@echo ""

$(SHOULD_COPY):
	@echo "Constructing $(SHOULD_COPY) script"
	@echo "import sys, os"                                              > $@
	@echo "filename = sys.argv[1]"                                     >> $@
	@echo "if os.path.abspath('$(srcdir)') == os.getcwd():"            >> $@
	@echo "    print 'false'; sys.exit()"                              >> $@
	@echo "if not os.path.exists(filename):"                           >> $@
	@echo "    print 'true'; sys.exit()"                               >> $@
	@echo "if (os.path.getmtime(filename) <"                           >> $@
	@echo "    os.path.getmtime(os.path.join('$(srcdir)',filename))):" >> $@
	@echo "    print 'true'; sys.exit()"                               >> $@
	@echo "print 'false'"                                              >> $@

all-local: $(SHOULD_COPY) $(EXAMPLE_SUITE)
	@if [ `$(PYTHON) $(SHOULD_COPY) $(SETPATH)` == 'true' ]; then \
	  echo "copying $(SETPATH)"; /bin/cp $(srcdir)/$(SETPATH) . ; fi
	@if [ `$(PYTHON) $(SHOULD_COPY) $(GRE__115)` == 'true' ]; then \
	  echo "copying $(GRE__115)"; /bin/cp $(srcdir)/$(GRE__115) . ; fi

ex%: $(SHOULD_COPY)
	@if [ `$(PYTHON) $(SHOULD_COPY) $@.py` == 'true' ]; then \
	  echo "copying $@.py"; /bin/cp $(srcdir)/$@.py . ; fi

run-examples: all-local
	@for ex_script in $(EXAMPLE_PYSCRIPTS); do if [ $$ex_script != '.py' ]; then \
	  echo; echo "Running $$ex_script"; echo; $$ex_script -t; fi; done

CLEANFILES = $(EXAMPLE_PYSCRIPTS) $(SHOULD_COPY) $(SETPATH) $(GRE__115) \
	     *.pyc *.dat *.mm *.xml
