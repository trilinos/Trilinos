# @HEADER
# ************************************************************************
# 
#              PyTrilinos: Python Interface to Trilinos
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
# Questions? Contact Bill Spotz (wfspotz@sandia.gov) 
# 
# ************************************************************************
# @HEADER

EXTRA_DIST = setup.py

# Configuration macros
@HAVE_MPI_TRUE@HAVE_MPI = 1
@HAVE_MPI_FALSE@HAVE_MPI =

@HAVE_PYTRILINOS_TEUCHOS_TRUE@ENABLE_TEUCHOS=true
@HAVE_PYTRILINOS_TEUCHOS_FALSE@ENABLE_TEUCHOS=

@HAVE_PYTRILINOS_THYRA_TRUE@ENABLE_THYRA=true
@HAVE_PYTRILINOS_THYRA_FALSE@ENABLE_THYRA=

@HAVE_PYTRILINOS_EPETRA_TRUE@ENABLE_EPETRA=true
@HAVE_PYTRILINOS_EPETRA_FALSE@ENABLE_EPETRA=

@HAVE_PYTRILINOS_TRIUTILS_TRUE@ENABLE_TRIUTILS=true
@HAVE_PYTRILINOS_TRIUTILS_FALSE@ENABLE_TRIUTILS=

@HAVE_PYTRILINOS_EPETRAEXT_TRUE@ENABLE_EPETRAEXT=true
@HAVE_PYTRILINOS_EPETRAEXT_FALSE@ENABLE_EPETRAEXT=

@HAVE_PYTRILINOS_AZTECOO_TRUE@ENABLE_AZTECOO=true
@HAVE_PYTRILINOS_AZTECOO_FALSE@ENABLE_AZTECOO=
@HAVE_PYTRILINOS_AZTECOO_TEUCHOS_TRUE@ENABLE_AZTECOO_TEUCHOS=true
@HAVE_PYTRILINOS_AZTECOO_TEUCHOS_FALSE@ENABLE_AZTECOO_TEUCHOS=

@HAVE_PYTRILINOS_GALERI_TRUE@ENABLE_GALERI=true
@HAVE_PYTRILINOS_GALERI_FALSE@ENABLE_GALERI=

@HAVE_PYTRILINOS_AMESOS_TRUE@ENABLE_AMESOS=true
@HAVE_PYTRILINOS_AMESOS_FALSE@ENABLE_AMESOS=

@HAVE_PYTRILINOS_IFPACK_TRUE@ENABLE_IFPACK=true
@HAVE_PYTRILINOS_IFPACK_FALSE@ENABLE_IFPACK=

@HAVE_PYTRILINOS_KOMPLEX_TRUE@ENABLE_KOMPLEX=true
@HAVE_PYTRILINOS_KOMPLEX_FALSE@ENABLE_KOMPLEX=

@HAVE_PYTRILINOS_ANASAZI_TRUE@ENABLE_ANASAZI=true
@HAVE_PYTRILINOS_ANASAZI_FALSE@ENABLE_ANASAZI=

@HAVE_PYTRILINOS_ML_TRUE@ENABLE_ML=true
@HAVE_PYTRILINOS_ML_FALSE@ENABLE_ML=

@HAVE_PYTRILINOS_NOX_TRUE@ENABLE_NOX=true
@HAVE_PYTRILINOS_NOX_FALSE@ENABLE_NOX=
@HAVE_PYTRILINOS_NOX_EPETRA_TRUE@ENABLE_NOX_EPETRA=true
@HAVE_PYTRILINOS_NOX_EPETRA_FALSE@ENABLE_NOX_EPETRA=

@HAVE_PYTRILINOS_LOCA_TRUE@ENABLE_LOCA=true
@HAVE_PYTRILINOS_LOCA_FALSE@ENABLE_LOCA=

# Default make target
all-local: package

# The following line helps the test harness recover from build errors.
all-am:
	@echo ""
	@echo "Trilinos package PyTrilinos subdirectory shared built successfully."
	@echo ""

package:
	$(PYTHON) $(srcdir)/setup.py build

install:
	$(PYTHON) $(srcdir)/setup.py install

uninstall:
	$(PYTHON) $(srcdir)/setup.py uninstall

clean-local:
	$(PYTHON) $(srcdir)/setup.py clean

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------
MAINTAINERCLEANFILES = Makefile.in 
