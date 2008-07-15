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

@BUILD_PYTRILINOS_TRUE@ENABLE_PYTRILINOS=true
@BUILD_PYTRILINOS_FALSE@ENABLE_PYTRILINOS=

@HAVE_PYTRILINOS_TEUCHOS_TRUE@ENABLE_TEUCHOS=true
@HAVE_PYTRILINOS_TEUCHOS_FALSE@ENABLE_TEUCHOS=

@HAVE_PYTRILINOS_RTOP_TRUE@ENABLE_RTOP=true
@HAVE_PYTRILINOS_RTOP_FALSE@ENABLE_RTOP=

@HAVE_PYTRILINOS_KOKKOS_TRUE@ENABLE_KOKKOS=true
@HAVE_PYTRILINOS_KOKKOS_FALSE@ENABLE_KOKKOS=

@HAVE_PYTRILINOS_EPETRA_TRUE@ENABLE_EPETRA=true
@HAVE_PYTRILINOS_EPETRA_FALSE@ENABLE_EPETRA=

@HAVE_PYTRILINOS_TRIUTILS_TRUE@ENABLE_TRIUTILS=true
@HAVE_PYTRILINOS_TRIUTILS_FALSE@ENABLE_TRIUTILS=

@HAVE_PYTRILINOS_TPETRA_TRUE@ENABLE_TPETRA=true
@HAVE_PYTRILINOS_TPETRA_FALSE@ENABLE_TPETRA=

@HAVE_PYTRILINOS_EPETRAEXT_TRUE@ENABLE_EPETRAEXT=true
@HAVE_PYTRILINOS_EPETRAEXT_FALSE@ENABLE_EPETRAEXT=

@HAVE_PYTRILINOS_THYRA_TRUE@ENABLE_THYRA=true
@HAVE_PYTRILINOS_THYRA_FALSE@ENABLE_THYRA=

@HAVE_PYTRILINOS_ISORROPIA_TRUE@ENABLE_ISORROPIA=true
@HAVE_PYTRILINOS_ISORROPIA_FALSE@ENABLE_ISORROPIA=

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

@HAVE_PYTRILINOS_CLAPS_TRUE@ENABLE_CLAPS=true
@HAVE_PYTRILINOS_CLAPS_FALSE@ENABLE_CLAPS=

@HAVE_PYTRILINOS_BELOS_TRUE@ENABLE_BELOS=true
@HAVE_PYTRILINOS_BELOS_FALSE@ENABLE_BELOS=

@HAVE_PYTRILINOS_ANASAZI_TRUE@ENABLE_ANASAZI=true
@HAVE_PYTRILINOS_ANASAZI_FALSE@ENABLE_ANASAZI=

@HAVE_PYTRILINOS_ML_TRUE@ENABLE_ML=true
@HAVE_PYTRILINOS_ML_FALSE@ENABLE_ML=

@HAVE_PYTRILINOS_MOERTEL_TRUE@ENABLE_MOERTEL=true
@HAVE_PYTRILINOS_MOERTEL_FALSE@ENABLE_MOERTEL=

@HAVE_PYTRILINOS_STRATIMIKOS_TRUE@ENABLE_STRATIMIKOS=true
@HAVE_PYTRILINOS_STRATIMIKOS_FALSE@ENABLE_STRATIMIKOS=

@HAVE_PYTRILINOS_MEROS_TRUE@ENABLE_MEROS=true
@HAVE_PYTRILINOS_MEROS_FALSE@ENABLE_MEROS=

@HAVE_PYTRILINOS_SACADO_TRUE@ENABLE_SACADO=true
@HAVE_PYTRILINOS_SACADO_FALSE@ENABLE_SACADO=

@HAVE_PYTRILINOS_NOX_TRUE@ENABLE_NOX=true
@HAVE_PYTRILINOS_NOX_FALSE@ENABLE_NOX=
@HAVE_PYTRILINOS_NOX_EPETRA_TRUE@ENABLE_NOX_EPETRA=true
@HAVE_PYTRILINOS_NOX_EPETRA_FALSE@ENABLE_NOX_EPETRA=

@HAVE_PYTRILINOS_LOCA_TRUE@ENABLE_LOCA=true
@HAVE_PYTRILINOS_LOCA_FALSE@ENABLE_LOCA=

@HAVE_PYTRILINOS_RYTHMOS_TRUE@ENABLE_RYTHMOS=true
@HAVE_PYTRILINOS_RYTHMOS_FALSE@ENABLE_RYTHMOS=

@HAVE_PYTRILINOS_MOOCHO_TRUE@ENABLE_MOOCHO=true
@HAVE_PYTRILINOS_MOOCHO_FALSE@ENABLE_MOOCHO=

# Teuchos
if HAVE_PYTRILINOS_TEUCHOS
include $(top_builddir)/../teuchos/Makefile.export.teuchos
endif

# RTOp
if HAVE_PYTRILINOS_RTOP
include $(top_builddir)/../rtop/Makefile.export.rtop
endif

# Epetra
if HAVE_PYTRILINOS_EPETRA
include $(top_builddir)/../epetra/Makefile.export.epetra
endif

# TriUtils
if HAVE_PYTRILINOS_TRIUTILS
include $(top_builddir)/../triutils/Makefile.export.triutils
endif

# Tpetra
if HAVE_PYTRILINOS_TPETRA
include $(top_builddir)/../tpetra/Makefile.export.tpetra
endif

# EpetraExt
if HAVE_PYTRILINOS_EPETRAEXT
include $(top_builddir)/../epetraext/Makefile.export.epetraext
endif

# Thyra
if HAVE_PYTRILINOS_THYRA
include $(top_builddir)/../thyra/Makefile.export.thyra
endif

# Isorropia
if HAVE_PYTRILINOS_ISORROPIA
include $(top_builddir)/../isorropia/Makefile.export.isorropia
endif

# AztecOO
if HAVE_PYTRILINOS_AZTECOO
include $(top_builddir)/../aztecoo/Makefile.export.aztecoo
endif

# Galeri
if HAVE_PYTRILINOS_GALERI
include $(top_builddir)/../galeri/Makefile.export.galeri
endif

# Amesos
if HAVE_PYTRILINOS_AMESOS
include $(top_builddir)/../amesos/Makefile.export.amesos
endif

# IFPACK
if HAVE_PYTRILINOS_IFPACK
include $(top_builddir)/../ifpack/Makefile.export.ifpack
endif

# Komplex
if HAVE_PYTRILINOS_KOMPLEX
include $(top_builddir)/../komplex/Makefile.export.komplex
endif

# Claps
if HAVE_PYTRILINOS_CLAPS
include $(top_builddir)/../claps/Makefile.export.claps
endif

# Belos
if HAVE_PYTRILINOS_BELOS
include $(top_builddir)/../belos/Makefile.export.belos
endif

# Anasazi
if HAVE_PYTRILINOS_ANASAZI
include $(top_builddir)/../anasazi/Makefile.export.anasazi
endif

# ML
if HAVE_PYTRILINOS_ML
include $(top_builddir)/../ml/Makefile.export.ml
endif

# Moertel
if HAVE_PYTRILINOS_MOERTEL
include $(top_builddir)/../moertel/Makefile.export.moertel
endif

# Stratimikos
if HAVE_PYTRILINOS_STRATIMIKOS
include $(top_builddir)/../stratimikos/Makefile.export.stratimikos
endif

# Meros
if HAVE_PYTRILINOS_MEROS
include $(top_builddir)/../meros/Makefile.export.meros
endif

# Sacado
if HAVE_PYTRILINOS_SACADO
include $(top_builddir)/../sacado/Makefile.export.sacado
endif

# NOX
if HAVE_PYTRILINOS_NOX
include $(top_builddir)/../nox/Makefile.export.nox
endif

# Rythmos
if HAVE_PYTRILINOS_RYTHMOS
include $(top_builddir)/../rythmos/Makefile.export.rythmos
endif

# Moocho
if HAVE_PYTRILINOS_MOOCHO
include $(top_builddir)/../moocho/Makefile.export.moocho
endif

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
	$(PYTHON) $(srcdir)/setup.py install --destdir=$(DESTDIR)

uninstall:
	$(PYTHON) $(srcdir)/setup.py uninstall

clean-local:
	$(PYTHON) $(srcdir)/setup.py clean

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------
MAINTAINERCLEANFILES = Makefile.in 
