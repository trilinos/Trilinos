# @HEADER
# ************************************************************************
# 
#       Ifpack: Object-Oriented Algebraic Preconditioner Package
#                 Copyright (2002) Sandia Corporation
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

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package ifpack subdirectory src built successfully."
	@echo ""

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside Ifpack. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#Ifpack_CrsGraph.h
#Ifpack_CrsRick.h
#
#
#

CORE_H = \
	$(srcdir)/Ifpack_ConfigDefs.h

CORE = 

UTIL_H = \
	$(srcdir)/Ifpack_Version.h

SERIAL_H = 

SERIAL = 

if HAVE_MPI

MPI_H = 

MPI = 

else

MPI_H =
MPI =

endif

#Ifpack does not yet support threads
#if HAVE_THREADS

MPISMP_H = 

MPISMP = 

#else

#MPISMP_H =
#MPISMP =

#endif

USERBASE_H = 

LANG_WRAPS_H = 

LANG_WRAPS = 

DIST_H = 
DIST = 

COMP_H = \
	$(srcdir)/Ifpack_IlukGraph.h           \
	$(srcdir)/Ifpack_CrsRiluk.h            \
	$(srcdir)/Ifpack_ILU.h                 \
	$(srcdir)/Ifpack_IKLU.h                \
	$(srcdir)/Ifpack_IKLU_Utils.h          \
	$(srcdir)/Ifpack_ILUT.h                \
	$(srcdir)/Ifpack_OverlapGraph.h        \
	$(srcdir)/Ifpack_ScalingType.h         \
	$(srcdir)/Ifpack_CrsIct.h              \
	$(srcdir)/Ifpack_ICT.h                 \
	$(srcdir)/Ifpack_IC.h                  \
	$(srcdir)/Ifpack_IC_Utils.h            \
	$(srcdir)/icrout_cholesky_mex.h        \
        $(srcdir)/ifp_parameters.h	       \
	$(srcdir)/Ifpack_Preconditioner.h      \
	$(srcdir)/Ifpack_CondestType.h         \
	$(srcdir)/Ifpack_Condest.h             \
	$(srcdir)/Ifpack.h                     \
	$(srcdir)/Ifpack_DiagPreconditioner.h  \
	$(srcdir)/Ifpack_PointRelaxation.h     \
	$(srcdir)/Ifpack_BlockRelaxation.h     \
	$(srcdir)/Ifpack_Container.h           \
	$(srcdir)/Ifpack_SparseContainer.h     \
	$(srcdir)/Ifpack_DenseContainer.h      \
	$(srcdir)/Ifpack_AdditiveSchwarz.h     \
	$(srcdir)/Ifpack_LocalFilter.h         \
	$(srcdir)/Ifpack_NodeFilter.h         \
	$(srcdir)/Ifpack_DropFilter.h          \
	$(srcdir)/Ifpack_SparsityFilter.h      \
	$(srcdir)/Ifpack_ReorderFilter.h       \
	$(srcdir)/Ifpack_SingletonFilter.h     \
	$(srcdir)/Ifpack_DiagonalFilter.h      \
	$(srcdir)/Ifpack_Graph.h               \
        $(srcdir)/Ifpack_Graph_Epetra_RowMatrix.h \
	$(srcdir)/Ifpack_Graph_Epetra_CrsGraph.h  \
	$(srcdir)/Ifpack_Partitioner.h            \
	$(srcdir)/Ifpack_OverlappingPartitioner.h \
	$(srcdir)/Ifpack_LinearPartitioner.h      \
	$(srcdir)/Ifpack_GreedyPartitioner.h      \
	$(srcdir)/Ifpack_METISPartitioner.h       \
	$(srcdir)/Ifpack_EquationPartitioner.h    \
	$(srcdir)/Ifpack_UserPartitioner.h        \
	$(srcdir)/Ifpack_Utils.h                  \
	$(srcdir)/Ifpack_HashTable.h              \
	$(srcdir)/Ifpack_Reordering.h             \
	$(srcdir)/Ifpack_RCMReordering.h          \
	$(srcdir)/Ifpack_SPARSKIT.h               \
	$(srcdir)/Ifpack_Chebyshev.h              \
	$(srcdir)/Ifpack_METISReordering.h        \
	$(srcdir)/Ifpack_OverlappingRowMatrix.h   \
	$(srcdir)/Ifpack_ValidParameters.h


COMP = \
	$(srcdir)/Ifpack_IlukGraph.cpp              \
	$(srcdir)/Ifpack_CrsRiluk.cpp               \
	$(srcdir)/Ifpack_ILU.cpp                    \
	$(srcdir)/Ifpack_IKLU.cpp                   \
	$(srcdir)/Ifpack_IKLU_Utils.cpp             \
	$(srcdir)/Ifpack_ILUT.cpp                   \
	$(srcdir)/Ifpack_CrsIct.cpp                 \
	$(srcdir)/Ifpack_ICT.cpp                    \
	$(srcdir)/Ifpack_IC.cpp                     \
	$(srcdir)/Ifpack_IC_Utils.cpp               \
	$(srcdir)/Ifpack_OverlapGraph.cpp           \
	$(srcdir)/icrout_cholesky_mex.c             \
	$(srcdir)/icrout_quicksort.c                \
	$(srcdir)/ifp_parameters.cpp                \
	$(srcdir)/Ifpack_DiagPreconditioner.cpp     \
	$(srcdir)/Ifpack_PointRelaxation.cpp        \
	$(srcdir)/Ifpack_Condest.cpp                \
	$(srcdir)/Ifpack.cpp                        \
	$(srcdir)/Ifpack_DenseContainer.cpp         \
	$(srcdir)/Ifpack_LocalFilter.cpp            \
	$(srcdir)/Ifpack_NodeFilter.cpp            \
	$(srcdir)/Ifpack_DropFilter.cpp             \
	$(srcdir)/Ifpack_SparsityFilter.cpp         \
	$(srcdir)/Ifpack_ReorderFilter.cpp          \
	$(srcdir)/Ifpack_SingletonFilter.cpp        \
	$(srcdir)/Ifpack_DiagonalFilter.cpp         \
	$(srcdir)/Ifpack_Graph_Epetra_RowMatrix.cpp \
	$(srcdir)/Ifpack_Graph_Epetra_CrsGraph.cpp  \
	$(srcdir)/Ifpack_Utils.cpp                  \
	$(srcdir)/Ifpack_Chebyshev.cpp              \
	$(srcdir)/Ifpack_OverlappingPartitioner.cpp \
	$(srcdir)/Ifpack_LinearPartitioner.cpp      \
	$(srcdir)/Ifpack_GreedyPartitioner.cpp      \
	$(srcdir)/Ifpack_METISPartitioner.cpp       \
	$(srcdir)/Ifpack_EquationPartitioner.cpp    \
	$(srcdir)/Ifpack_UserPartitioner.cpp        \
	$(srcdir)/Ifpack_RCMReordering.cpp          \
	$(srcdir)/Ifpack_METISReordering.cpp        \
	$(srcdir)/Ifpack_OverlappingRowMatrix.cpp   \
	$(srcdir)/Ifpack_ValidParameters.cpp


# Conditionally compile in Amesos support files

EXTRA_AMESOS_COMP_H = \
	$(srcdir)/Ifpack_AMDReordering.h \
        $(srcdir)/Ifpack_Amesos.h

EXTRA_AMESOS_COMP = \
	$(srcdir)/Ifpack_AMDReordering.cpp \
        $(srcdir)/Ifpack_Amesos.cpp

if HAVE_IFPACK_AMESOS
AMESOS_COMP_H = $(EXTRA_AMESOS_COMP_H)
AMESOS_COMP = $(EXTRA_AMESOS_COMP)
else
AMESOS_COMP_H = 
AMESOS_COMP = 
endif

# ------------------------------------------------------------------------
# Ifpack library specifications
# ------------------------------------------------------------------------

IFPACK_LIB = libifpack.a


IFPACK_H = \
	$(CORE_H) \
	$(UTIL_H) \
	$(SERIAL_H) \
	$(USERBASE_H) \
	$(LANG_WRAPS_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(AMESOS_COMP_H) \
	$(MPI_H) \
	$(MPISMP_H)

libifpack_a_SOURCES = \
	$(CORE) \
	$(SERIAL) \
	$(LANG_WRAPS) \
	$(DIST) \
	$(COMP) \
	$(AMESOS_COMP) \
	$(MPI) \
	$(MPISMP)

EXTRA_libifpack_a_SOURCES = \
	$(EXTRA_DIST) \
	$(EXTRA_DIST_H) 

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libifpack_a_AR = $(ALTERNATE_AR)
else

libifpack_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Includes for required and optional library dependencies
# NOTE:  Amesos path is conditionally added for AMD
# ------------------------------------------------------------------------
include $(top_builddir)/Makefile.export.ifpack

if HAVE_IFPACK_AMESOS
IFPACK_INCLUDES_EXTRA = \
	$(IFPACK_INCLUDES) \
	-I$(top_srcdir)/../amesos/src/SuiteSparse/AMD/Include \
	-I$(top_srcdir)/../amesos/src/SuiteSparse/UFconfig
else
IFPACK_INCLUDES_EXTRA = $(IFPACK_INCLUDES)
endif

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(IFPACK_INCLUDES_EXTRA))
else
EXPORT_INCLUDES = $(IFPACK_INCLUDES_EXTRA)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libifpack_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(IFPACK_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(IFPACK_H)

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Ifpack_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(includedir)
nodist_execinclude_HEADERS = Ifpack_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 
