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

EXTRA_DIST = Amesos.i Anasazi.i AztecOO.i \
	     Epetra.i Epetra_Base.i Epetra_Comm.i Epetra_Graphs.i \
	     Epetra_Maps.i Epetra_Operators.i Epetra_SerialDense.i Epetra_Vectors.i \
	     EpetraExt.i Galeri.i IFPACK.i \
	     LOCA.__init__.i LOCA.Abstract.i LOCA.Bifurcation.i LOCA.Chan.i \
	     LOCA.Continuation.i LOCA.Epetra.i LOCA.Homotopy.i LOCA.LAPACK.i \
	     LOCA.Multicontinuation.i LOCA.TimeDependent.i \
	     ML.i \
	     NOX.__init__.i NOX.Abstract.i NOX.Epetra.__init__.i \
	     NOX.Epetra.Interface.i NOX.Solver.i NOX.StatusTest.i \
	     Teuchos.i Teuchos_ParameterList_ext.i \
	     Thyra.i TriuTils.i \
	     numpy.i setup.py SharedUtils.py

@HAVE_MPI_TRUE@HAVE_MPI = 1
@HAVE_MPI_FALSE@HAVE_MPI =

# Default make target
all-local: build

# Include swig rules
include $(top_srcdir)/config/swig.am

# For the package-enabling macros that follow, I would have preferred
# all info (interface file names, python proxy file names, wrapper code
# file names and wrapper header file names) to derive from a single
# source (interface file name, for example).  This would be possible
# with GNU make substitution rules, but this is not portable.  So I
# require a list of interface files and a list of proxy files for each
# package.

# Teuchos
if HAVE_PYTRILINOS_TEUCHOS
include $(top_builddir)/../teuchos/Makefile.export.teuchos
TEUCHOS_INTERFACES = Teuchos.i
TEUCHOS_PROXIES = PyTrilinos/Teuchos.py
endif

# Thyra
if HAVE_PYTRILINOS_THYRA
include $(top_builddir)/../thyra/Makefile.export.thyra
THYRA_INTERFACES = Thyra.i
THYRA_PROXIES = PyTrilinos/Thyra.py
endif

# Epetra
if HAVE_PYTRILINOS_EPETRA
include $(top_builddir)/../epetra/Makefile.export.epetra
EPETRA_INTERFACES = Epetra.i
EPETRA_PROXIES = PyTrilinos/Epetra.py
endif

# TriUtils
if HAVE_PYTRILINOS_TRIUTILS
include $(top_builddir)/../triutils/Makefile.export.triutils
TRIUTILS_INTERFACES = TriUtils.i
TRIUTILS_PROXIES = PyTrilinos/TriUtils.py
endif

# EpetraExt
if HAVE_PYTRILINOS_EPETRAEXT
include $(top_builddir)/../epetraext/Makefile.export.epetraext
EPETRAEXT_INTERFACES = EpetraExt.i
EPETRAEXT_PROXIES = PyTrilinos/EpetraExt.py
endif

# AztecOO
if HAVE_PYTRILINOS_AZTECOO
include $(top_builddir)/../aztecoo/Makefile.export.aztecoo
AZTECOO_INTERFACES = AztecOO.i
AZTECOO_PROXIES = PyTrilinos/AztecOO.py
endif

# Galeri
if HAVE_PYTRILINOS_GALERI
include $(top_builddir)/../galeri/Makefile.export.galeri
GALERI_INTERFACES = Galeri.i
GALERI_PROXIES = PyTrilinos/Galeri.py
endif

# Amesos
if HAVE_PYTRILINOS_AMESOS
include $(top_builddir)/../amesos/Makefile.export.amesos
AMESOS_INTERFACES = Amesos.i
AMESOS_PROXIES = PyTrilinos/Amesos.py
endif

# IFPACK
if HAVE_PYTRILINOS_IFPACK
include $(top_builddir)/../ifpack/Makefile.export.ifpack
IFPACK_INTERFACES = IFPACK.i
IFPACK_PROXIES = PyTrilinos/IFPACK.py
endif

# Anasazi
if HAVE_PYTRILINOS_ANASAZI
include $(top_builddir)/../anasazi/Makefile.export.anasazi
ANASAZI_INTERFACES = Anasazi.i
ANASAZI_PROXIES = PyTrilinos/Anasazi.py
endif

# ML
if HAVE_PYTRILINOS_ML
include $(top_builddir)/../ml/Makefile.export.ml
ML_INTERFACES = ML.i
ML_PROXIES = PyTrilinos/ML.py
endif

# NOX
if HAVE_PYTRILINOS_NOX
include $(top_builddir)/../nox/Makefile.export.nox
NOX_BASE_INTERFACES = NOX.__init__.i NOX.Abstract.i NOX.Solver.i NOX.StatusTest.i
NOX_BASE_PROXIES = PyTrilinos/NOX/__init__.py PyTrilinos/NOX/Abstract.py \
                   PyTrilinos/NOX/Solver.py   PyTrilinos/NOX/StatusTest.py
if HAVE_PYTRILINOS_NOX_EPETRA
NOX_EPETRA_INTERFACES = NOX.Epetra.__init__.i NOX.Epetra.Interface.i
NOX_EPETRA_PROXIES = PyTrilinos/NOX/Epetra/__init__.py PyTrilinos/NOX/Epetra/Interface.py
endif
NOX_INTERFACES = $(NOX_BASE_INTERFACES) $(NOX_EPETRA_INTERFACES) $(NOX_LAPACK_INTERFACES)
NOX_PROXIES = $(NOX_BASE_PROXIES) $(NOX_EPETRA_PROXIES) $(NOX_LAPACK_PROXIES)
endif

# LOCA
if HAVE_PYTRILINOS_NOX
if HAVE_PYTRILINOS_LOCA
LOCA_BASE_INTERFACES =  LOCA.__init__.i LOCA.Abstract.i LOCA.Bifurcation.i LOCA.Chan.i \
	                LOCA.Continuation.i LOCA.Homotopy.i LOCA.MultiContinuation.i \
	                LOCA.TimeDependent.i
LOCA_BASE_PROXIES =  PyTrilinos/LOCA/__init__.py PyTrilinos/LOCA/Abstract.py \
	             PyTrilinos/LOCA/Bifurcation.py PyTrilinos/LOCA/Chan.py \
	             PyTrilinos/LOCA/Continuation.py PyTrilinos/LOCA/Homotopy.py \
	             PyTrilinos/LOCA/MultiContinuation.py PyTrilinos/LOCA/TimeDependent.py
if HAVE_PYTRILINOS_NOX_EPETRA
LOCA_EPETRA_INTERFACES = LOCA.Epetra.i
LOCA_EPETRA_PROXIES = PyTrilinos/LOCA/Epetra.py
endif
#if HAVE_PYTRILINOS_NOX_LAPACK
#LOCA_LAPACK_INTERFACES = LOCA.LAPACK.i
#LOCA_LAPACK_PROXIES = PyTrilinos/LOCA/LAPACK.py
#endif
LOCA_INTERFACES = $(LOCA_BASE_INTERFACES) $(LOCA_EPETRA_INTERFACES) $(LOCA_LAPACK_INTERFACES)
LOCA_PROXIES = $(LOCA_BASE_PROXIES) $(LOCA_EPETRA_PROXIES) $(LOCA_LAPACK_PROXIES)
endif
endif

# PyTrilinos
INTERFACES = $(TEUCHOS_INTERFACES)  $(THYRA_INTERFACES)     $(EPETRA_INTERFACES)  \
	     $(TRIUTILS_INTERFACES) $(EPETRAEXT_INTERFACES) $(AZTECOO_INTERFACES) \
	     $(GALERI_INTERFACES)   $(AMESOS_INTERFACES)    $(IFPACK_INTERFACES)  \
	     $(ANASAZI_INTERFACES)  $(ML_INTERFACES)        $(NOX_INTERFACES)     \
	     $(LOCA_INTERFACES)
PROXIES  =   $(TEUCHOS_PROXIES)  $(THYRA_PROXIES)     $(EPETRA_PROXIES)  \
	     $(TRIUTILS_PROXIES) $(EPETRAEXT_PROXIES) $(AZTECOO_PROXIES) \
	     $(GALERI_PROXIES)   $(AMESOS_PROXIES)    $(IFPACK_PROXIES)  \
	     $(ANASAZI_PROXIES)  $(ML_PROXIES)        $(NOX_PROXIES)     \
	     $(LOCA_INTERFACES)
WRAPPERS   = $(INTERFACES:.i=_wrap.cpp)
DIR_HDRS   = $(INTERFACES:.i=_wrap.h)

# The following line helps the test harness recover from build errors.
all-am:
	@echo ""
	@echo "Trilinos package PyTrilinos subdirectory src built successfully."
	@echo ""

build: $(WRAPPERS)
	$(PYTHON) $(srcdir)/setup.py build

install: build
	$(PYTHON) $(srcdir)/setup.py install

uninstall:
	$(PYTHON) $(srcdir)/setup.py uninstall

clean-local:
	$(PYTHON) $(srcdir)/setup.py clean

CLEANFILES += $(PROXIES) $(WRAPPERS) $(DIR_HDRS)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------
MAINTAINERCLEANFILES = Makefile.in 
