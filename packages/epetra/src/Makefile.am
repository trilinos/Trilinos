# @HEADER
# ************************************************************************
# 
#               Epetra: Linear Algebra Services Package 
#                 Copyright (2001) Sandia Corporation
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


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside Epetra. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

CORE_H = \
	$(srcdir)/Epetra_BLAS.h \
	$(srcdir)/Epetra_LAPACK.h \
	$(srcdir)/Epetra_BLAS_wrappers.h \
	$(srcdir)/Epetra_LAPACK_wrappers.h \
	$(srcdir)/Epetra_Comm.h \
	$(srcdir)/Epetra_Object.h \
	$(srcdir)/Epetra_ConfigDefs.h


CORE = \
	$(srcdir)/Epetra_BLAS.cpp \
	$(srcdir)/Epetra_LAPACK.cpp \
	$(srcdir)/Epetra_Object.cpp
	
UTIL_H = \
	$(srcdir)/Epetra_Flops.h \
	$(srcdir)/Epetra_Time.h \
	$(srcdir)/Epetra_Util.h \
	$(srcdir)/Epetra_DataAccess.h \
	$(srcdir)/Epetra_HashTable.h \
	$(srcdir)/Epetra_RowMatrixTransposer.h \
	$(srcdir)/Epetra_Data.h \
	$(srcdir)/Epetra_Version.h

UTIL = \
	$(srcdir)/Epetra_Flops.cpp \
	$(srcdir)/Epetra_Time.cpp \
	$(srcdir)/Epetra_Util.cpp \
	$(srcdir)/Epetra_RowMatrixTransposer.cpp \
	$(srcdir)/Epetra_Data.cpp

ENUM_H = \
	$(srcdir)/Epetra_CombineMode.h \
	$(srcdir)/Epetra_DataAccess.h

UNIPROC_H = \
	$(srcdir)/Epetra_SerialDenseOperator.h  \
	$(srcdir)/Epetra_SerialDenseMatrix.h  \
	$(srcdir)/Epetra_SerialSpdDenseSolver.h \
	$(srcdir)/Epetra_SerialDenseSolver.h   \
	$(srcdir)/Epetra_SerialSymDenseMatrix.h \
	$(srcdir)/Epetra_SerialDenseVector.h \
	$(srcdir)/Epetra_SerialDenseSVD.h \
	$(srcdir)/Epetra_IntSerialDenseMatrix.h \
	$(srcdir)/Epetra_IntSerialDenseVector.h

UNIPROC = \
	$(srcdir)/Epetra_SerialDenseMatrix.cpp  \
	$(srcdir)/Epetra_SerialSpdDenseSolver.cpp \
	$(srcdir)/Epetra_SerialDenseSolver.cpp   \
	$(srcdir)/Epetra_SerialSymDenseMatrix.cpp \
	$(srcdir)/Epetra_SerialDenseVector.cpp \
	$(srcdir)/Epetra_SerialDenseSVD.cpp \
	$(srcdir)/Epetra_IntSerialDenseMatrix.cpp  \
	$(srcdir)/Epetra_IntSerialDenseVector.cpp 

SERIAL_H = \
	$(srcdir)/Epetra_SerialComm.h  \
	$(srcdir)/Epetra_SerialCommData.h  \
	$(srcdir)/Epetra_SerialDistributor.h

SERIAL = \
	$(srcdir)/Epetra_SerialComm.cpp  \
	$(srcdir)/Epetra_SerialCommData.cpp  \
	$(srcdir)/Epetra_SerialDistributor.cpp

EXTRA_MPI_H = \
	$(srcdir)/Epetra_MpiComm.h  \
	$(srcdir)/Epetra_MpiCommData.h  \
	$(srcdir)/Epetra_MpiDistributor.h 

EXTRA_MPI = \
	$(srcdir)/Epetra_MpiComm.cpp  \
	$(srcdir)/Epetra_MpiCommData.cpp  \
	$(srcdir)/Epetra_MpiDistributor.cpp

if HAVE_MPI

MPI_H = \
	$(EXTRA_MPI_H)
MPI = \
	$(EXTRA_MPI)

else

MPI_H =
MPI =

endif

EXTRA_MPISMP_H = \
	$(srcdir)/Epetra_MpiComm.h \
	$(srcdir)/Epetra_MpiCommData.h \
	$(srcdir)/Epetra_MpiSmpComm.h \
	$(srcdir)/Epetra_MpiSmpCommData.h

EXTRA_MPISMP = \
	$(srcdir)/Epetra_MpiComm.cpp \
	$(srcdir)/Epetra_MpiCommData.cpp \
	$(srcdir)/Epetra_MpiSmpComm.cpp \
	$(srcdir)/Epetra_MpiSmpCommData.cpp

if HAVE_THREADS

MPISMP_H = \
	$(EXTRA_MPISMP_H)
MPISMP = \
	$(EXTRA_MPISMP)

else

MPISMP_H =
MPISMP =

endif



USERBASE_H = \
	$(srcdir)/Epetra_Operator.h \
	$(srcdir)/Epetra_RowMatrix.h


LANG_WRAPS_H = \
	$(srcdir)/Epetra_C_wrappers.h \
	$(srcdir)/Epetra_Fortran_wrappers.h

LANG_WRAPS = \
	$(srcdir)/Epetra_C_wrappers.cpp \
	$(srcdir)/Epetra_Fortran_wrappers.cpp

DIST_H = \
	$(srcdir)/Epetra_BasicDirectory.h \
	$(srcdir)/Epetra_BlockMap.h \
	$(srcdir)/Epetra_BlockMapData.h \
	$(srcdir)/Epetra_Directory.h \
	$(srcdir)/Epetra_SrcDistObject.h \
	$(srcdir)/Epetra_DistObject.h \
	$(srcdir)/Epetra_Distributor.h \
	$(srcdir)/Epetra_Export.h \
	$(srcdir)/Epetra_Import.h \
	$(srcdir)/Epetra_LocalMap.h \
	$(srcdir)/Epetra_Map.h \
	$(srcdir)/Epetra_MapColoring.h \
        $(srcdir)/Epetra_OffsetIndex.h

DIST = \
	$(srcdir)/Epetra_BasicDirectory.cpp \
	$(srcdir)/Epetra_BlockMap.cpp \
	$(srcdir)/Epetra_BlockMapData.cpp \
	$(srcdir)/Epetra_DistObject.cpp \
	$(srcdir)/Epetra_Export.cpp \
	$(srcdir)/Epetra_Import.cpp \
	$(srcdir)/Epetra_LocalMap.cpp \
	$(srcdir)/Epetra_Map.cpp \
	$(srcdir)/Epetra_MapColoring.cpp \
        $(srcdir)/Epetra_OffsetIndex.cpp


COMP_H = \
	$(srcdir)/Epetra_CompObject.h \
	$(srcdir)/Epetra_CrsGraph.h \
	$(srcdir)/Epetra_CrsGraphData.h \
	$(srcdir)/Epetra_CrsMatrix.h \
	$(srcdir)/Epetra_LinearProblem.h \
	$(srcdir)/Epetra_MultiVector.h \
	$(srcdir)/Epetra_VbrMatrix.h \
	$(srcdir)/Epetra_Vector.h \
	$(srcdir)/Epetra_CrsSingletonFilter.h \
	$(srcdir)/Epetra_FECrsMatrix.h \
	$(srcdir)/Epetra_FEVbrMatrix.h \
	$(srcdir)/Epetra_FEVector.h \
	$(srcdir)/Epetra_IntVector.h

COMP = \
	$(srcdir)/Epetra_CompObject.cpp \
	$(srcdir)/Epetra_CrsGraph.cpp \
	$(srcdir)/Epetra_CrsGraphData.cpp \
	$(srcdir)/Epetra_CrsMatrix.cpp \
	$(srcdir)/Epetra_LinearProblem.cpp \
	$(srcdir)/Epetra_MultiVector.cpp \
	$(srcdir)/Epetra_VbrMatrix.cpp \
	$(srcdir)/Epetra_Vector.cpp \
	$(srcdir)/Epetra_CrsSingletonFilter.cpp \
	$(srcdir)/Epetra_FECrsMatrix.cpp \
	$(srcdir)/Epetra_FEVbrMatrix.cpp \
	$(srcdir)/Epetra_FEVector.cpp \
	$(srcdir)/Epetra_IntVector.cpp


# ------------------------------------------------------------------------
# Epetra library specifications
# ------------------------------------------------------------------------

EPETRA_LIB = libepetra.a


EPETRA_H = \
	$(CORE_H) \
	$(UTIL_H) \
	$(ENUM_H) \
	$(UNIPROC_H) \
	$(SERIAL_H) \
	$(USERBASE_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H) \
	$(MPISMP_H)

libepetra_a_SOURCES = \
	$(CORE) \
	$(UTIL) \
	$(UNIPROC) \
	$(SERIAL) \
	$(DIST) \
	$(COMP) \
	$(MPI) \
	$(MPISMP)

EXTRA_libepetra_a_SOURCES = \
	$(EXTRA_MPI) \
	$(EXTRA_MPI_H) \
	$(EXTRA_MPISMP) \
	$(EXTRA_MPISMP_H)

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libepetra_a_AR = $(ALTERNATE_AR)
else

libepetra_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# On alpha-dec machines, the extra .o's are needed for Epetra only.
# ------------------------------------------------------------------------

libepetra_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(EPETRA_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(EPETRA_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Epetra_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = Epetra_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






