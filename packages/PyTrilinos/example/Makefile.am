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

EXTRA_DIST = setpath.src gre__115.mx \
	     exAmesos_Factory.src exAmesos_Solvers.src \
	     exAnasazi_BlockDavidson.src exAztecOO.src exAztecOO_Operator.src \
	     exAztecOO_RowMatrix.src exAztecOO_BasicRowMatrix.src exEpetra.src \
	     exEpetra_Comm.src exEpetra_CrsMatrix_Easy.src \
	     exEpetra_CrsMatrix_Efficient.src exEpetra_FECrsMatrix_Easy.src \
	     exEpetra_ImportExport.src exEpetraExt_IO_MatrixMarket.src \
	     exEpetraExt_IO_XML.src exEpetraExt_MatrixMatrix.src exGaleri.src \
	     exGaleri_ReadHB.src exIFPACK.src exKomplex_Simple.src exMLAPI.src \
	     exMLAPI_Simple.src exMLAPI_AztecOO.src exMLAPI_Iterate.src \
	     exMLAPI_PyMatrix.src exMLAPI_Smoother.src exNOX_1Dfem.src \
	     exTeuchos.src exThyra.src exLOCA_Chan.src

MAINTAINERCLEANFILES = Makefile.in

## #######################################################################
## Define the suite of examples
## #######################################################################

# Teuchos examples
if HAVE_PYTRILINOS_TEUCHOS
TEUCHOS_EXAMPLES = exTeuchos.py
endif

# # Thyra examples
# if HAVE_PYTRILINOS_THYRA
# if HAVE_PYTRILINOS_EPETRA
# THYRA_EXAMPLES = exThyra.py
# endif
# endif

# Epetra examples
if HAVE_PYTRILINOS_EPETRA
EPETRA_EXAMPLES = exEpetra.py exEpetra_Comm.py exEpetra_CrsMatrix_Easy.py \
	          exEpetra_CrsMatrix_Efficient.py exEpetra_FECrsMatrix_Easy.py \
	          exEpetra_ImportExport.py
endif

# TriUtils examples
if HAVE_PYTRILINOS_TRIUTILS
TRIUTILS_EXAMPLES = 
endif

# EpetraExt examples
if HAVE_PYTRILINOS_EPETRAEXT
if HAVE_PYTRILINOS_EPETRA
EPETRAEXT_EXAMPLES = exEpetraExt_IO_MatrixMarket.py exEpetraExt_IO_XML.py \
	             exEpetraExt_MatrixMatrix.py
endif
endif

# AztecOO examples
if HAVE_PYTRILINOS_AZTECOO
if HAVE_PYTRILINOS_EPETRA
if HAVE_PYTRILINOS_GALERI
AZTECOO_EX1 = exAztecOO.py
endif
AZTECOO_EX2 = exAztecOO_Operator.py exAztecOO_RowMatrix.py \
	      exAztecOO_BasicRowMatrix.py
AZTECOO_EXAMPLES = $(AZTECOO_EX1) $(AZTECOO_EX2)
endif
endif

# Galeri examples
if HAVE_PYTRILINOS_GALERI
if HAVE_PYTRILINOS_EPETRA
GALERI_EXAMPLES = exGaleri.py exGaleri_ReadHB.py
endif
endif

# Amesos examples
if HAVE_PYTRILINOS_AMESOS
if HAVE_PYTRILINOS_EPETRA
AMESOS_EXAMPLES = exAmesos_Factory.py exAmesos_Solvers.py
endif
endif

# IFPACK examples
if HAVE_PYTRILINOS_IFPACK
if HAVE_PYTRILINOS_AZTECOO
if HAVE_PYTRILINOS_GALERI
if HAVE_PYTRILINOS_TRIUTILS
if HAVE_PYTRILINOS_EPETRA
IFPACK_EXAMPLES = exIFPACK.py
endif
endif
endif
endif
endif

# Komplex examples
if HAVE_PYTRILINOS_KOMPLEX
if HAVE_PYTRILINOS_AZTECOO
if HAVE_PYTRILINOS_EPETRA
KOMPLEX_EXAMPLES = exKomplex_Simple.py
endif
endif
endif

# Anasazi examples
if HAVE_PYTRILINOS_ANASAZI
if HAVE_PYTRILINOS_GALERI
ANASAZI_EXAMPLES = exAnasazi_BlockDavidson.py
endif
endif

# ML examples
if HAVE_PYTRILINOS_ML
if HAVE_PYTRILINOS_EPETRA
if HAVE_MPI
ML_EX1 = exMLAPI_PyMatrix.py exMLAPI_Smoother.py
else
ML_EX1 = exMLAPI.py exMLAPI_Iterate.py exMLAPI_PyMatrix.py exMLAPI_Smoother.py
endif
if HAVE_PYTRILINOS_AZTECOO
if HAVE_MPI
ML_EX2 = 
else
ML_EX2 = exMLAPI_AztecOO.py
endif
if HAVE_PYTRILINOS_GALERI
ML_EX3 = exMLAPI_Simple.py
endif
endif
ML_EXAMPLES = $(ML_EX1) $(ML_EX2) $(ML_EX3)
endif
endif

# NOX examples
if HAVE_PYTRILINOS_NOX
if HAVE_PYTRILINOS_NOX_EPETRA
if HAVE_MPI
NOX_EXAMPLES = 
else
NOX_EXAMPLES = exNOX_1Dfem.py
endif
endif
endif

# LOCA examples
if HAVE_PYTRILINOS_NOX
if HAVE_PYTRILINOS_LOCA
if HAVE_PYTRILINOS_NOX_EPETRA
LOCA_EXAMPLES = exLOCA_Chan.py
endif
endif
endif

# Collect the package examples into a suite of examples
EXAMPLE_SUITE = $(TEUCHOS_EXAMPLES)  $(THYRA_EXAMPLES)     $(EPETRA_EXAMPLES) \
	        $(TRIUTILS_EXAMPLES) $(EPETRAEXT_EXAMPLES) $(AZTECOO_EXAMPLES) \
	        $(GALERI_EXAMPLES)   $(AMESOS_EXAMPLES)    $(IFPACK_EXAMPLES) \
	        $(KOMPLEX_EXAMPLES)  $(ANASAZI_EXAMPLES)   $(ML_EXAMPLES) \
	        $(NOX_EXAMPLES)      $(LOCA_EXAMPLES)

SETPATH  = setpath.py
GRE__115 = gre__115.rua

# The following line helps the test harness recover from build errors.
all:

all-local: $(SETPATH) $(EXAMPLE_SUITE) $(GRE__115)

%.py: %.src
	@echo "copying $@"; /bin/cp $< $@

%.rua: %.mx
	@echo "copying $@"; /bin/cp $< $@

run-examples: all-local
	@for ex_script in $(EXAMPLE_SUITE); do echo; echo "Running $$ex_script"; \
	 echo; $$ex_script -t; done

CLEANFILES = $(EXAMPLE_SUITE) $(SETPATH) $(GRE__115) *.pyc *.dat *.mm *.xml
