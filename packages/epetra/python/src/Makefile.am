# @HEADER
# ************************************************************************
# 
#              PyTrilinos.Epetra: Python Interface to Epetra
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

# Include export Makefiles
include $(top_builddir)/Makefile.export.epetra
include $(top_builddir)/../PyTrilinos/Makefile.export.pytrilinos

# Include swig rules
include $(top_srcdir)/config/swig.am
@HAVE_MPI_TRUE@MPI_FLAG = -DHAVE_MPI

SWIG_CXXFLAGS = $(EPETRA_INCLUDES)     \
	        $(PYTRILINOS_INCLUDES) \
	        -I$(srcdir)            \
		$(MPI_FLAG)

EXTRA_DIST = Epetra.i                             \
	     Epetra_Base.i                        \
	     Epetra_Comm.i                        \
	     Epetra_Graphs.i                      \
	     Epetra_Maps.i                        \
	     Epetra_Operators.i                   \
	     Epetra_SerialDense.i                 \
	     Epetra_Vectors.i                     \
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
	     Epetra_PyOperator.h                  \
	     Epetra_PyRowMatrix.h                 \
	     setup.py

# Default make target
all-am: build

# Currently, this temporary directory is required to be present for the
# setup.py script to work correctly
tempdir = $(shell $(PYTHON) -c "import sys,os; from distutils.util import get_platform; print os.path.join('build','temp.%s-%s'%(get_platform(),sys.version[0:3]))")

$(tempdir):
	@echo "Creating directory $@"
	mkdir -p $@

# Define the SWIG interface files and the resulting wrapper files
INTERFACES = Epetra.i
WRAPPERS   = $(INTERFACES:.i=_wrap.cpp)
PROXIES    = $(INTERFACES:.i=.py      )

build: $(WRAPPERS) $(tempdir)
	$(PYTHON) $(srcdir)/setup.py build

install: build
	$(PYTHON) $(srcdir)/setup.py install --prefix=$(PYTHON_PREFIX)

clean-local:
	$(PYTHON) $(srcdir)/setup.py clean
	$(RM) -r ../../packages

CLEANFILES += $(WRAPPERS) $(PROXIES)
