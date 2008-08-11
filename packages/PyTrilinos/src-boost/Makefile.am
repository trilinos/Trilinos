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

EXTRA_DIST = MakefileVariables.py pkg_info.py PyTrilinosExtension.py setup.py \
	     SharedUtils.py \
	     Teuchos_wrap.cpp Teuchos.py

# Configuration macros
@HAVE_MPI_TRUE@HAVE_MPI = 1
@HAVE_MPI_FALSE@HAVE_MPI =

@HAVE_PYTRILINOS_TEUCHOS_TRUE@ENABLE_TEUCHOS=true
@HAVE_PYTRILINOS_TEUCHOS_FALSE@ENABLE_TEUCHOS=

#@HAVE_PYTRILINOS_THYRA_TRUE@ENABLE_THYRA=true
#@HAVE_PYTRILINOS_THYRA_FALSE@ENABLE_THYRA=
ENABLE_THYRA=

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

@HAVE_PYTRILINOS_ANASAZI_TRUE@ENABLE_ANASAZI=true
@HAVE_PYTRILINOS_ANASAZI_FALSE@ENABLE_ANASAZI=

@HAVE_PYTRILINOS_ML_TRUE@ENABLE_ML=true
@HAVE_PYTRILINOS_ML_FALSE@ENABLE_ML=

@HAVE_PYTRILINOS_NOX_TRUE@ENABLE_NOX=true
@HAVE_PYTRILINOS_NOX_FALSE@ENABLE_NOX=
@HAVE_PYTRILINOS_NOX_EPETRA_TRUE@ENABLE_NOX_EPETRA=true
@HAVE_PYTRILINOS_NOX_EPETRA_FALSE@ENABLE_NOX_EPETRA=

#@HAVE_PYTRILINOS_LOCA_TRUE@ENABLE_LOCA=true
#@HAVE_PYTRILINOS_LOCA_FALSE@ENABLE_LOCA=
ENABLE_LOCA=

# Default make target
all-local: package

# Currently, this temporary directory is required to be present for the
# setup.py script to work correctly
tempdir = $(shell $(PYTHON) -c "import sys,os; from distutils.util import get_platform; print os.path.join('build','temp.%s-%s'%(get_platform(),sys.version[0:3]))")
libdir = $(shell $(PYTHON) -c "import sys,os; from distutils.util import get_platform; print os.path.join('build','lib.%s-%s'%(get_platform(),sys.version[0:3]),'PyTrilinos')")

$(tempdir):
	@echo "creating directory $@"
	@mkdir -p $@

# For the package-enabling macros that follow, I would have preferred
# all info (wrapper file names and proxy file names) to derive from a
# single source (interface file name, for example).  This would be
# possible with GNU make substitution rules, but this is not portable.
# So I require a list of wrapper files, a list of proxy file names and
# a list of module names for each package.

# Teuchos
if HAVE_PYTRILINOS_TEUCHOS
include $(top_builddir)/../teuchos/Makefile.export.teuchos
TEUCHOS_WRAPPERS = Teuchos_wrap.cpp
TEUCHOS_PROXIES  = PyTrilinos/Teuchos.py
endif

# # Thyra
# if HAVE_PYTRILINOS_THYRA
# include $(top_builddir)/../thyra/Makefile.export.thyra
# endif

# Epetra
if HAVE_PYTRILINOS_EPETRA
include $(top_builddir)/../epetra/Makefile.export.epetra
endif

# TriUtils
if HAVE_PYTRILINOS_TRIUTILS
include $(top_builddir)/../triutils/Makefile.export.triutils
endif

# EpetraExt
if HAVE_PYTRILINOS_EPETRAEXT
include $(top_builddir)/../epetraext/Makefile.export.epetraext
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

# Anasazi
if HAVE_PYTRILINOS_ANASAZI
include $(top_builddir)/../anasazi/Makefile.export.anasazi
endif

# ML
if HAVE_PYTRILINOS_ML
include $(top_builddir)/../ml/Makefile.export.ml
endif

# NOX
if HAVE_PYTRILINOS_NOX
include $(top_builddir)/../nox/Makefile.export.nox
endif

# # LOCA
# if HAVE_PYTRILINOS_NOX
# if HAVE_PYTRILINOS_LOCA
# endif
# endif

# PyTrilinos wrappers and proxies
WRAPPERS = $(TEUCHOS_WRAPPERS)
PROXIES  = $(TEUCHOS_PROXIES)

# The following line helps the test harness recover from build errors.
all-am:
	@echo ""
	@echo "Trilinos package PyTrilinos subdirectory src built successfully."
	@echo ""

PyTrilinos/%.py: $(srcdir)/%.py
	@echo "copying $< to $@"
	cp $< $@

package: $(tempdir) $(PROXIES)
	$(PYTHON) $(srcdir)/setup.py build

install:
	$(PYTHON) $(srcdir)/setup.py install

uninstall:
	$(PYTHON) $(srcdir)/setup.py uninstall

clean-local:
	$(RM) -r $(libdir)
	$(RM) -r $(tempdir)/$(srcdir)
	$(PYTHON) $(srcdir)/setup.py clean

CLEANFILES = $(PROXIES)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------
MAINTAINERCLEANFILES = Makefile.in 
