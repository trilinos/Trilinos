# @HEADER
# ************************************************************************
# 
#     EpetraExt: Epetra Extended - Linear Algebra Services Package 
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
# might someday be needed by some other code outside EpetraExt. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

CORE_H = \
     $(srcdir)/EpetraExt_ConfigDefs.h \
     $(srcdir)/EpetraExt_MatrixMatrix.h \
     $(srcdir)/EpetraExt_Version.h

CORE = \
     $(srcdir)/EpetraExt_MatrixMatrix.cpp

EXTRA_TRANSFORM_H = \
     $(srcdir)/transform/EpetraExt_Transform.h \
     $(srcdir)/transform/EpetraExt_Transform_Composite.h \
     $(srcdir)/transform/EpetraExt_Permutation.h \
     $(srcdir)/transform/EpetraExt_LPTrans_From_GraphTrans.h \
     $(srcdir)/transform/EpetraExt_LPTrans_From_MatrixTrans.h \
     $(srcdir)/transform/EpetraExt_View_CrsGraph.h \
     $(srcdir)/transform/EpetraExt_View_CrsMatrix.h \
     $(srcdir)/transform/EpetraExt_View_MultiVector.h \
     $(srcdir)/transform/EpetraExt_Transpose_CrsGraph.h \
     $(srcdir)/transform/EpetraExt_Transpose_RowMatrix.h \
     $(srcdir)/transform/EpetraExt_Scale_LinearProblem.h \
     $(srcdir)/transform/EpetraExt_CrsSingletonFilter_LinearProblem.h \
     $(srcdir)/transform/EpetraExt_Overlap_CrsGraph.h \
     $(srcdir)/transform/EpetraExt_Reindex_CrsMatrix.h \
     $(srcdir)/transform/EpetraExt_Reindex_MultiVector.h \
     $(srcdir)/transform/EpetraExt_Reindex_LinearProblem.h \
     $(srcdir)/transform/EpetraExt_SolverMap_CrsMatrix.h \
     $(srcdir)/transform/EpetraExt_SolverMap_LinearProblem.h
     
EXTRA_TRANSFORM = \
     $(srcdir)/transform/EpetraExt_Permutation.cpp \
     $(srcdir)/transform/EpetraExt_LPTrans_From_GraphTrans.cpp \
     $(srcdir)/transform/EpetraExt_LPTrans_From_MatrixTrans.cpp \
     $(srcdir)/transform/EpetraExt_View_CrsGraph.cpp \
     $(srcdir)/transform/EpetraExt_View_CrsMatrix.cpp \
     $(srcdir)/transform/EpetraExt_View_MultiVector.cpp \
     $(srcdir)/transform/EpetraExt_Transpose_CrsGraph.cpp \
     $(srcdir)/transform/EpetraExt_Transpose_RowMatrix.cpp \
     $(srcdir)/transform/EpetraExt_Scale_LinearProblem.cpp \
     $(srcdir)/transform/EpetraExt_CrsSingletonFilter_LinearProblem.cpp \
     $(srcdir)/transform/EpetraExt_Overlap_CrsGraph.cpp \
     $(srcdir)/transform/EpetraExt_Reindex_CrsMatrix.cpp \
     $(srcdir)/transform/EpetraExt_Reindex_MultiVector.cpp \
     $(srcdir)/transform/EpetraExt_Reindex_LinearProblem.cpp \
     $(srcdir)/transform/EpetraExt_SolverMap_CrsMatrix.cpp \
     $(srcdir)/transform/EpetraExt_SolverMap_LinearProblem.cpp

EXTRA_TRANSFORM_EXPERIMENTAL_H = \
     $(srcdir)/transform/EpetraExt_SymmRCM_CrsGraph.h \
     $(srcdir)/transform/EpetraExt_StaticCondensation_LinearProblem.h \
     $(srcdir)/transform/EpetraExt_AMD_CrsGraph.h

EXTRA_TRANSFORM_EXPERIMENTAL = \
     $(srcdir)/transform/EpetraExt_SymmRCM_CrsGraph.cpp \
     $(srcdir)/transform/EpetraExt_StaticCondensation_LinearProblem.cpp \
     $(srcdir)/transform/EpetraExt_AMD_CrsGraph.cpp

if BUILD_TRANSFORM

if BUILD_EXPERIMENTAL

TRANSFORM_H = \
     $(EXTRA_TRANSFORM_H) $(EXTRA_TRANSFORM_EXPERIMENTAL_H)

TRANSFORM = \
     $(EXTRA_TRANSFORM) $(EXTRA_TRANSFORM_EXPERIMENTAL)

else

TRANSFORM_H = \
     $(EXTRA_TRANSFORM_H)

TRANSFORM = \
     $(EXTRA_TRANSFORM)

endif

else

TRANSFORM_H =

TRANSFORM =

endif
     
EXTRA_BTF_H = \
     $(srcdir)/btf/EpetraExt_BTF_CrsGraph.h \
     $(srcdir)/btf/EpetraExt_BTF_CrsMatrix.h \
     $(srcdir)/btf/EpetraExt_BTF_LinearProblem.h \
     $(srcdir)/btf/EpetraExt_BlockJacobi_LinearProblem.h

EXTRA_BTF = \
     $(srcdir)/btf/EpetraExt_BTF_CrsGraph.cpp \
     $(srcdir)/btf/EpetraExt_BTF_CrsMatrix.cpp \
     $(srcdir)/btf/EpetraExt_BTF_LinearProblem.cpp \
     $(srcdir)/btf/EpetraExt_BlockJacobi_LinearProblem.cpp \
     $(srcdir)/btf/pothen/concmp.f \
     $(srcdir)/btf/pothen/genbtf.f \
     $(srcdir)/btf/pothen/corsum.f \
     $(srcdir)/btf/pothen/finsum.f \
     $(srcdir)/btf/pothen/fnrsum.f \
     $(srcdir)/btf/pothen/izero.f \
     $(srcdir)/btf/pothen/maxmatch.f \
     $(srcdir)/btf/pothen/rectblk.f \
     $(srcdir)/btf/pothen/mmc13e.f \
     $(srcdir)/btf/pothen/prtivs.f \
     $(srcdir)/btf/pothen/strchk.f \
     $(srcdir)/btf/pothen/mattrans.f

if BUILD_BTF

BTF_H = \
     $(EXTRA_BTF_H)

BTF = \
     $(EXTRA_BTF)

else

BTF_H =

BTF =

endif

EXTRA_INOUT_H = \
     $(srcdir)/inout/EpetraExt_mmio.h \
     $(srcdir)/inout/EpetraExt_RowMatrixOut.h \
     $(srcdir)/inout/EpetraExt_MultiVectorOut.h \
     $(srcdir)/inout/EpetraExt_VectorOut.h \
     $(srcdir)/inout/EpetraExt_BlockMapOut.h \
     $(srcdir)/inout/EpetraExt_BlockMapIn.h \
     $(srcdir)/inout/EpetraExt_CrsMatrixIn.h \
     $(srcdir)/inout/EpetraExt_MultiVectorIn.h \
     $(srcdir)/inout/EpetraExt_VectorIn.h

EXTRA_INOUT = \
     $(srcdir)/inout/EpetraExt_mmio.cpp \
     $(srcdir)/inout/EpetraExt_RowMatrixOut.cpp \
     $(srcdir)/inout/EpetraExt_MultiVectorOut.cpp \
     $(srcdir)/inout/EpetraExt_VectorOut.cpp \
     $(srcdir)/inout/EpetraExt_BlockMapOut.cpp \
     $(srcdir)/inout/EpetraExt_BlockMapIn.cpp \
     $(srcdir)/inout/EpetraExt_CrsMatrixIn.cpp \
     $(srcdir)/inout/EpetraExt_MultiVectorIn.cpp \
     $(srcdir)/inout/EpetraExt_VectorIn.cpp

if BUILD_INOUT

INOUT_H = \
    $(EXTRA_INOUT_H)

INOUT = \
    $(EXTRA_INOUT)

else

INOUT_H =

INOUT =

endif

EXTRA_ZOLTAN_H = \
     $(srcdir)/zoltan/EpetraExt_Zoltan_CrsGraph.h \
     $(srcdir)/zoltan/EpetraExt_ZoltanOrder_CrsGraph.h \
     $(srcdir)/zoltan/EpetraExt_ZoltanQuery.h \
     $(srcdir)/zoltan/Zoltan_LoadBalance.h \
     $(srcdir)/zoltan/Zoltan_MigrationContainer.h \
     $(srcdir)/zoltan/Zoltan_MigrationFunctions.h \
     $(srcdir)/zoltan/Zoltan_MigrationObject.h \
     $(srcdir)/zoltan/Zoltan_QueryContainer.h \
     $(srcdir)/zoltan/Zoltan_QueryFunctions.h \
     $(srcdir)/zoltan/Zoltan_QueryObject.h

EXTRA_ZOLTAN = \
     $(srcdir)/zoltan/EpetraExt_Zoltan_CrsGraph.cpp \
     $(srcdir)/zoltan/EpetraExt_ZoltanOrder_CrsGraph.cpp \
     $(srcdir)/zoltan/EpetraExt_ZoltanQuery.cpp \
     $(srcdir)/zoltan/Zoltan_LoadBalance.cpp \
     $(srcdir)/zoltan/Zoltan_MigrationContainer.cpp \
     $(srcdir)/zoltan/Zoltan_MigrationFunctions.cpp \
     $(srcdir)/zoltan/Zoltan_MigrationObject.cpp \
     $(srcdir)/zoltan/Zoltan_QueryContainer.cpp \
     $(srcdir)/zoltan/Zoltan_QueryFunctions.cpp \
     $(srcdir)/zoltan/Zoltan_QueryObject.cpp

#EXTRA_ZOLTAN_EXPERIMENTAL_H = \
#     $(srcdir)/zoltan/EpetraExt_ZoltanMpiComm.h \
#     $(srcdir)/zoltan/EpetraExt_ZoltanMpiCommData.h \
#     $(srcdir)/zoltan/EpetraExt_ZoltanMpiDistributor.h
EXTRA_ZOLTAN_EXPERIMENTAL_H = 
#EXTRA_ZOLTAN_EXPERIMENTAL = \
#     $(srcdir)/zoltan/EpetraExt_ZoltanMpiComm.cpp \
#     $(srcdir)/zoltan/EpetraExt_ZoltanMpiCommData.cpp \
#     $(srcdir)/zoltan/EpetraExt_ZoltanMpiDistributor.cpp
EXTRA_ZOLTAN_EXPERIMENTAL = 

if BUILD_ZOLTAN

if BUILD_EXPERIMENTAL

ZOLTAN_H = \
     $(EXTRA_ZOLTAN_H) $(EXTRA_ZOLTAN_EXPERIMENTAL_H)

ZOLTAN = \
     $(EXTRA_ZOLTAN) $(EXTRA_ZOLTAN_EXPERIMENTAL)

else

ZOLTAN_H = \
     $(EXTRA_ZOLTAN_H)

ZOLTAN = \
     $(EXTRA_ZOLTAN)

endif

else

ZOLTAN_H =

ZOLTAN =

endif

EXTRA_COLORING_H = \
     $(srcdir)/coloring/EpetraExt_MapColoring.h \
     $(srcdir)/coloring/EpetraExt_MapColoringIndex.h

EXTRA_COLORING = \
     $(srcdir)/coloring/EpetraExt_MapColoring.cpp \
     $(srcdir)/coloring/EpetraExt_MapColoringIndex.cpp

if BUILD_COLORING

COLORING_H = \
     $(EXTRA_COLORING_H)

COLORING = \
     $(EXTRA_COLORING)

else

COLORING_H =

COLORING =

endif

EXTRA_MATLAB_H = \
    $(srcdir)/matlab/EpetraExt_MatlabEngine.h \
    $(srcdir)/matlab/EpetraExt_PutMultiVector.h \
	$(srcdir)/matlab/EpetraExt_PutRowMatrix.h

EXTRA_MATLAB = \
    $(srcdir)/matlab/EpetraExt_MatlabEngine.cpp \
    $(srcdir)/matlab/EpetraExt_PutMultiVector.cpp \
	$(srcdir)/matlab/EpetraExt_PutRowMatrix.cpp

if BUILD_MATLAB

MATLAB_H = \
   $(EXTRA_MATLAB_H)

MATLAB = \
   $(EXTRA_MATLAB)

else

MATLAB_H =

MATLAB =

endif

# ------------------------------------------------------------------------
# EpetraExt inout library specifications
# ------------------------------------------------------------------------

libepetraext_a_SOURCES = \
     $(CORE) \
     $(TRANSFORM) \
     $(BTF) \
     $(INOUT) \
     $(ZOLTAN) \
     $(COLORING) \
     $(MATLAB)

#If these are not listed, the files will not be included in
#the distribution tarball.
EXTRA_libepetraext_a_SOURCES = 

# libraries to be installed in $(libexecdir)
EPETRAEXT_LIB = libepetraext.a

EPETRAEXT_H = \
     $(CORE_H) \
     $(TRANSFORM_H) \
     $(BTF_H) \
     $(INOUT_H) \
     $(ZOLTAN_H) \
     $(COLORING_H) \
     $(MATLAB_H)

#We need the following to support separate build directories
AM_CPPFLAGS = \
     -I$(top_builddir)/../epetra/src \
     -I$(top_srcdir)/../epetra/src \
     -I$(top_builddir)/../epetraext/src \
     -I$(top_srcdir)/../epetraext/src \
     -I$(top_builddir)/../epetraext/src/transform \
     -I$(top_srcdir)/../epetraext/src/transform \
     -I$(top_builddir)/../epetraext/src/btf \
     -I$(top_srcdir)/../epetraext/src/btf \
     -I$(top_builddir)/../epetraext/src/inout \
     -I$(top_srcdir)/../epetraext/src/inout \
     -I$(top_builddir)/../epetraext/src/zoltan \
     -I$(top_srcdir)/../epetraext/src/zoltan \
     -I$(top_builddir)/../epetraext/src/coloring \
     -I$(top_srcdir)/../epetraext/src/coloring \
     -I$(top_builddir)/../epetraext/src/matlab \
     -I$(top_srcdir)/../epetraext/src/matlab
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libepetraext_a_AR = $(ALTERNATE_AR)

else

libepetraext_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# On alpha-dec machines, the extra .o's are needed for Epetra only.
# ------------------------------------------------------------------------

libepetraext_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(EPETRAEXT_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(EPETRAEXT_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Epetra_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = EpetraExt_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Files to be deleted by 'make clean'
# ------------------------------------------------------------------------

MOSTLYCLEANFILES = $(TRILINOS_TEMPLATE_OBJS)

