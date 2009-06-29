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

EXTRA_DIST = Amesos.i                             \
	     Anasazi.i                            \
	     AztecOO.i                            \
	     Epetra.i                             \
	     EpetraExt.i                          \
	     Epetra_Base.i                        \
	     Epetra_Comm.i                        \
	     Epetra_Graphs.i                      \
	     Epetra_Maps.i                        \
	     Epetra_NumPyFEVector.cpp             \
	     Epetra_NumPyFEVector.h               \
	     Epetra_NumPyIntSerialDenseMatrix.cpp \
	     Epetra_NumPyIntSerialDenseMatrix.h   \
	     Epetra_NumPyIntSerialDenseVector.cpp \
	     Epetra_NumPyIntSerialDenseVector.h   \
	     Epetra_NumPyIntVector.cpp            \
	     Epetra_NumPyIntVector.h              \
	     Epetra_NumPyMultiVector.cpp          \
	     Epetra_NumPyMultiVector.h            \
	     Epetra_NumPySerialDenseMatrix.cpp    \
	     Epetra_NumPySerialDenseMatrix.h      \
	     Epetra_NumPySerialDenseVector.cpp    \
	     Epetra_NumPySerialDenseVector.h      \
	     Epetra_NumPyVector.cpp               \
	     Epetra_NumPyVector.h                 \
	     Epetra_Operators.i                   \
	     Epetra_RowMatrix_Utils.i             \
	     Epetra_SerialDense.i                 \
	     Epetra_Vectors.i                     \
             FILEstream.cpp                       \
             FILEstream.h                         \
	     Galeri.i                             \
	     IFPACK.i                             \
	     Komplex.i                            \
	     LOCA.Abstract.i                      \
	     LOCA.Bifurcation.i                   \
	     LOCA.Chan.i                          \
	     LOCA.Continuation.i                  \
	     LOCA.Epetra.Interface.i              \
	     LOCA.Epetra.__init__.i               \
	     LOCA.Homotopy.i                      \
	     LOCA.LAPACK.i                        \
	     LOCA.MultiContinuation.i             \
	     LOCA.TurningPoint.i                  \
	     LOCA.Pitchfork.i                     \
	     LOCA.TimeDependent.i                 \
	     LOCA.Extended.i     \
	     LOCA.BorderedSystem.i   \
	     LOCA.__init__.i                      \
	     ML.i                                 \
	     MLAPI_PyMatrix.h                     \
	     NOX.Abstract.i                       \
	     NOX.Epetra.Interface.i               \
	     NOX.Epetra.__init__.i                \
	     NOX.Solver.i                         \
	     NOX.StatusTest.i                     \
	     NOX.__init__.i                       \
	     NumPyImporter.h                      \
	     Pliris.i                             \
	     PythonException.cpp                  \
	     PythonException.h                    \
	     Teuchos.i                            \
	     Teuchos_Traits.i                     \
	     Teuchos_Comm.i                       \
	     Teuchos_Epetra.i                     \
	     Teuchos_ParameterList.i              \
	     Teuchos_Traits.i                     \
	     Teuchos_XML.i                        \
	     Teuchos_PythonParameter.cpp          \
	     Teuchos_PythonParameter.h            \
	     Thyra.i                              \
	     TriUtils.i                           \
	     UserArray.patch                      \
	     numpy.i                              \
	     numpy_include.h                      \
	     pyfragments.swg.in                   \
	     setup.py                             \
	     scripts/ParamConvert.py

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
@HAVE_PYTRILINOS_EPETRAEXT_HDF5_TRUE@ENABLE_EPETRAEXT_HDF5=true
@HAVE_PYTRILINOS_EPETRAEXT_HDF5_FALSE@ENABLE_EPETRAEXT_HDF5=

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

@HAVE_PYTRILINOS_PLIRIS_TRUE@ENABLE_PLIRIS=true
@HAVE_PYTRILINOS_PLIRIS_FALSE@ENABLE_PLIRIS=

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

# Include swig rules
include $(top_srcdir)/config/swig.am
if HAVE_MPI
SWIG_OPTS += -DHAVE_MPI
endif

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
SWIG_OPTS += -DHAVE_TEUCHOS
endif

# # Thyra
# if HAVE_PYTRILINOS_THYRA
# include $(top_builddir)/../thyra/Makefile.export.thyra
# THYRA_INTERFACES = Thyra.i
# THYRA_PROXIES = PyTrilinos/Thyra.py
# endif

# Epetra
if HAVE_PYTRILINOS_EPETRA
include $(top_builddir)/../epetra/Makefile.export.epetra
EPETRA_INTERFACES = Epetra.i
EPETRA_PROXIES = PyTrilinos/Epetra.py
SWIG_OPTS += -DHAVE_EPETRA
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
if HAVE_PYTRILINOS_EPETRAEXT_HDF5
SWIG_OPTS += -DHAVE_EPETRAEXT_HDF5
endif
endif

# AztecOO
if HAVE_PYTRILINOS_AZTECOO
include $(top_builddir)/../aztecoo/Makefile.export.aztecoo
AZTECOO_INTERFACES = AztecOO.i
AZTECOO_PROXIES = PyTrilinos/AztecOO.py
# if HAVE_PYTRILINOS_AZTECOO_TEUCHOS
# SWIG_OPTS += -DHAVE_AZTECOO_TEUCHOS 
# endif
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

# Komplex
if HAVE_PYTRILINOS_KOMPLEX
include $(top_builddir)/../komplex/Makefile.export.komplex
KOMPLEX_INTERFACES = Komplex.i
KOMPLEX_PROXIES = PyTrilinos/Komplex.py
endif

# Anasazi
if HAVE_PYTRILINOS_ANASAZI
include $(top_builddir)/../anasazi/Makefile.export.anasazi
ANASAZI_INTERFACES = Anasazi.i
ANASAZI_PROXIES = PyTrilinos/Anasazi.py
endif

# Pliris
if HAVE_PYTRILINOS_PLIRIS
include $(top_builddir)/../pliris/Makefile.export.pliris
PLIRIS_INTERFACES = Pliris.i
PLIRIS_PROXIES = PyTrilinos/Pliris.py
endif

# ML
if HAVE_PYTRILINOS_ML
include $(top_builddir)/../ml/Makefile.export.ml
ML_INTERFACES = ML.i
ML_PROXIES = PyTrilinos/ML.py
SWIG_OPTS += -DMLAPI_LC
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
SWIG_OPTS += -DHAVE_NOX_EPETRA
endif

# LOCA
if HAVE_PYTRILINOS_NOX
if HAVE_PYTRILINOS_LOCA
LOCA_BASE_INTERFACES = LOCA.__init__.i LOCA.Abstract.i LOCA.Bifurcation.i \
	               LOCA.TurningPoint.i LOCA.Continuation.i            \
	               LOCA.Homotopy.i LOCA.MultiContinuation.i           \
	               LOCA.Pitchfork.i LOCA.TimeDependent.i LOCA.Hopf.i  \
	               LOCA.Extended.i LOCA.BorderedSystem.i
LOCA_BASE_PROXIES = PyTrilinos/LOCA/__init__.py          \
	            PyTrilinos/LOCA/Abstract.py          \
	            PyTrilinos/LOCA/TurningPoint.i       \
	            PyTrilinos/LOCA/Bifurcation.py       \
	            PyTrilinos/LOCA/Continuation.py      \
	            PyTrilinos/LOCA/Pitchfork.i          \
	            PyTrilinos/LOCA/Homotopy.py          \
	            PyTrilinos/LOCA/MultiContinuation.py \
	            PyTrilinos/LOCA/Hopf.i               \
	            PyTrilinos/LOCA/TimeDependent.py     \
	            PyTrilinos/LOCA/Extended.py          \
	            PyTrilinos/LOCA/BorderedSystem.py
if HAVE_PYTRILINOS_NOX_EPETRA
LOCA_EPETRA_INTERFACES = LOCA.Epetra.__init__.i LOCA.Epetra.Interface.i
LOCA_EPETRA_PROXIES = PyTrilinos/LOCA/Epetra/__init__.py PyTrilinos/LOCA/Epetra/Interface.py
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
	     $(KOMPLEX_INTERFACES)  $(ANASAZI_INTERFACES)   $(PLIRIS_INTERFACES)  \
	     $(ML_INTERFACES)       $(NOX_INTERFACES)       $(LOCA_INTERFACES)
PROXIES  =   $(TEUCHOS_PROXIES)  $(THYRA_PROXIES)     $(EPETRA_PROXIES)  \
	     $(TRIUTILS_PROXIES) $(EPETRAEXT_PROXIES) $(AZTECOO_PROXIES) \
	     $(GALERI_PROXIES)   $(AMESOS_PROXIES)    $(IFPACK_PROXIES)  \
	     $(KOMPLEX_PROXIES)  $(ANASAZI_PROXIES)   $(PLIRIS_PROXIES)  \
	     $(ML_PROXIES)       $(NOX_PROXIES)       $(LOCA_INTERFACES)
WRAPPERS   = $(INTERFACES:.i=_wrap.cpp)
DIR_HDRS   = $(INTERFACES:.i=_wrap.h)

# The following line helps the test harness recover from build errors.
all-am:
	@echo ""
	@echo "Trilinos package PyTrilinos subdirectory src built successfully."
	@echo ""

SUBDIRS = PyTrilinos

package: $(WRAPPERS)
	$(PYTHON) $(srcdir)/setup.py build

install: package
	if test -z '$(DESTDIR)'; then \
	$(PYTHON) $(srcdir)/setup.py install --prefix $(PYTHON_PREFIX); \
	else \
	$(PYTHON) $(srcdir)/setup.py install --root $(DESTDIR); \
	fi

uninstall:
	$(PYTHON) $(srcdir)/setup.py uninstall

clean-local:
	$(PYTHON) $(srcdir)/setup.py clean

CLEANFILES += $(PROXIES) $(WRAPPERS) $(DIR_HDRS)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------
MAINTAINERCLEANFILES = Makefile.in 
