# @HEADER
# ************************************************************************
# 
#               Thyra: Interfaces and Support Code for the Interoperability of Abstract Numerical Algorithms
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
# @HEADER


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
# There is, at present, no distinction between most of the categories, 
# other than those which start with EXTRA_.  All categories which do not
# start with EXTRA_ are included in all builds.  Categories which start
# with EXTRA_ are included in some builds (typically controlled by 
# --with-xyz or --without-xyz)
#

UTILS_H = \
  $(srcdir)/utils/Thyra_Range1D.hpp

UTILS = \
  $(srcdir)/utils/Thyra_Range1D.cpp

ANA_RTOP_BASE_H = \
  $(srcdir)/ana/rtop/base/RTOpPack_MPI_apply_op.hpp\
  $(srcdir)/ana/rtop/base/RTOpPack_RTOpServer.hpp\
  $(srcdir)/ana/rtop/base/RTOpPack_RTOpTDecl.hpp\
  $(srcdir)/ana/rtop/base/RTOpPack_TOpUnaryFuncPtrDecl.hpp\
  $(srcdir)/ana/rtop/base/RTOpPack_MPI_apply_op_decl.hpp\
  $(srcdir)/ana/rtop/base/RTOpPack_RTOpServerDecl.hpp\
  $(srcdir)/ana/rtop/base/RTOpPack_RTOpTHelpers.hpp\
  $(srcdir)/ana/rtop/base/RTOpPack_Types.hpp\
  $(srcdir)/ana/rtop/base/RTOpPack_RTOpT.hpp\
  $(srcdir)/ana/rtop/base/RTOpPack_TOpUnaryFuncPtr.hpp

ANA_RTOP_BASE = \
  $(srcdir)/ana/rtop/base/RTOpPack_MPI_apply_op.cpp

ANA_RTOP_BASE_C_H = \
  $(srcdir)/ana/rtop/base/RTOp_config.h\
  $(srcdir)/ana/rtop/base/RTOp_mpi.h\
  $(srcdir)/ana/rtop/base/RTOp_parallel_helpers.h

ANA_RTOP_BASE_C = \
  $(srcdir)/ana/rtop/base/RTOp_mpi.c\
  $(srcdir)/ana/rtop/base/RTOp_parallel_helpers.c

ANA_RTOP_OPS_LIB_H = \
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_ROpCountNanInf.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_ROpNorm1.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_TOpAssignVectors.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_ROpDotProd.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_ROpNorm2.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_TOpEleWiseDivide.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_ROpGetSubVector.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_ROpNormInf.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_TOpEleWiseProd.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_ROpMax.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_ROpSum.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_TOpLinearCombination.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_ROpMaxIndex.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_ROpWeightedNorm2.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_TOpRandomize.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_ROpMaxIndexLessThanBound.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_TOpAXPY.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_TOpReciprocal.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_ROpMin.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_TOpAbs.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_TOpScaleVector.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_ROpMinIndex.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_TOpAddScalar.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_TOpSetSubVector.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_ROpMinIndexGreaterThanBound.hpp\
  $(srcdir)/ana/rtop/ops_lib/RTOpPack_TOpAssignScalar.hpp

ANA_OPERATOR_VECTOR_INTERFACES_H = \
  $(srcdir)/ana/operator_vector/interfaces/Thyra_AssertOp.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_EuclideanLinearOpBase.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_EuclideanLinearOpBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_EuclideanScalarProd.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_EuclideanScalarProdDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_LinearOpBase.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_LinearOpBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_LinearOpScalarProd.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_LinearOpScalarProdDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_MultiVectorBase.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_MultiVectorBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_MultiVectorCols.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_MultiVectorColsDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_MultiVectorStdOps.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_MultiVectorStdOpsDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_OpBase.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_OpBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_ScalarProdBase.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_ScalarProdBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_ScalarProdVectorSpaceBase.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_ScalarProdVectorSpaceBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SerialMultiVectorBase.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SerialMultiVectorBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SerialMultiVectorStd.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SerialMultiVectorStdDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SerialVectorBase.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SerialVectorBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SerialVectorSpaceBase.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SerialVectorSpaceBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SerialVectorSpaceFactoryStd.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SerialVectorSpaceFactoryStdDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SerialVectorSpaceStd.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SerialVectorSpaceStdDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SerialVectorStd.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SerialVectorStdDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SingleRhsEuclideanLinearOpBase.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SingleRhsEuclideanLinearOpBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SingleRhsLinearOpBase.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SingleRhsLinearOpBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SingleScalarEuclideanLinearOpBase.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SingleScalarEuclideanLinearOpBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SingleScalarLinearOpBase.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_SingleScalarLinearOpBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_OperatorVectorTypes.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_VectorBase.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_VectorBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_VectorMultiVector.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_VectorMultiVectorDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_VectorSpaceBase.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_VectorSpaceBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_VectorSpaceFactoryBase.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_VectorSpaceFactoryBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_VectorStdOps.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_VectorStdOpsDecl.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_apply_op_helper.hpp\
  $(srcdir)/ana/operator_vector/interfaces/Thyra_apply_op_helper_decl.hpp

ANA_OPERATOR_VECTOR_EXTENDED_INTERFACES_H = \
  $(srcdir)/ana/operator_vector/extended_interfaces/Thyra_LinearOpWithSolveBase.hpp\
  $(srcdir)/ana/operator_vector/extended_interfaces/Thyra_LinearOpWithSolveBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/extended_interfaces/Thyra_LinearOpWithSolveFactoryBase.hpp\
  $(srcdir)/ana/operator_vector/extended_interfaces/Thyra_VectorSpaceConverterBase.hpp\
  $(srcdir)/ana/operator_vector/extended_interfaces/Thyra_ProductMultiVectorBase.hpp\
  $(srcdir)/ana/operator_vector/extended_interfaces/Thyra_ProductVectorBase.hpp\
  $(srcdir)/ana/operator_vector/extended_interfaces/Thyra_ProductVectorSpaceBase.hpp

ANA_OPERATOR_VECTOR_ANA_SUPPORT_H = \
  $(srcdir)/ana/operator_vector/ana_support/Thyra_DiagonalLinearOp.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_DiagonalLinearOpDecl.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_ExplicitMultiVectorView.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_ExplicitVectorView.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_LinearOpTester.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_LinearOpTesterDecl.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_MultiVectorAllocator.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_MultiVectorStdOpsTester.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_MultiVectorStdOpsTesterDecl.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_MultiVectorTester.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_MultiVectorTesterDecl.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_MultiplicativeLinearOp.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_MultiplicativeLinearOpDecl.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_ScaledAdjointLinearOp.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_ScaledAdjointLinearOpDecl.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_ProductVector.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_ProductVectorDecl.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_ProductVectorSpace.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_ProductVectorSpaceDecl.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_TestingTools.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_TestingToolsDecl.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_VectorSpaceTester.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_VectorSpaceTesterDecl.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_VectorStdOpsTester.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_VectorStdOpsTesterDecl.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_VectorTester.hpp\
  $(srcdir)/ana/operator_vector/ana_support/Thyra_VectorTesterDecl.hpp

ANA_OPERATOR_VECTOR_ADAPTER_SUPPORT_H = \
  $(srcdir)/ana/operator_vector/adapter_support/Thyra_SingleRhsLinearOpWithSolveBase.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/Thyra_SingleRhsLinearOpWithSolveBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/Thyra_SingleScalarLinearOpWithSolveBase.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/Thyra_SingleScalarLinearOpWithSolveBaseDecl.hpp

ANA_OPERATOR_VECTOR_ADAPTER_SUPPORT_SERIAL_H = \
  $(srcdir)/ana/operator_vector/adapter_support/serial/Thyra_SerialVectorSpaceConverterBase.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/serial/Thyra_SerialVectorSpaceConverterBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/serial/Thyra_SerialVectorSpaceConverterStd.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/serial/Thyra_SerialVectorSpaceConverterStdDecl.hpp\
  $(srcdir)/ana/operator_vector/extended_interfaces/Thyra_VectorSpaceConverterBase.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/serial/Thyra_SerialLinearOpBase.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/serial/Thyra_SerialLinearOpBaseDecl.hpp

ANA_OPERATOR_VECTOR_ADAPTER_SUPPORT_MPI_H = \
  $(srcdir)/ana/operator_vector/adapter_support/mpi/Thyra_MPILinearOpBase.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/mpi/Thyra_MPILinearOpBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/mpi/Thyra_MPIMultiVectorBase.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/mpi/Thyra_MPIMultiVectorBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/mpi/Thyra_MPIMultiVectorStd.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/mpi/Thyra_MPIMultiVectorStdDecl.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/mpi/Thyra_MPIVectorBase.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/mpi/Thyra_MPIVectorBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/mpi/Thyra_MPIVectorSpaceBase.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/mpi/Thyra_MPIVectorSpaceBaseDecl.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/mpi/Thyra_MPIVectorSpaceFactoryStd.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/mpi/Thyra_MPIVectorSpaceFactoryStdDecl.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/mpi/Thyra_MPIVectorSpaceStd.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/mpi/Thyra_MPIVectorSpaceStdDecl.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/mpi/Thyra_MPIVectorStd.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/mpi/Thyra_MPIVectorStdDecl.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/mpi/Thyra_MultiVectorSerialization.hpp\
  $(srcdir)/ana/operator_vector/adapter_support/mpi/Thyra_MultiVectorSerializationDecl.hpp

#Specify locations of include files within Trilinos

include $(top_builddir)/Makefile.export.thyra

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(THYRA_INCLUDES))
else
EXPORT_INCLUDES = $(THYRA_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# --------------------------------------------------------------------
# Thyra library specifications
# --------------------------------------------------------------------

THYRA_LIB = libthyra.a

THYRA_H = \
  $(srcdir)/Thyra_ConfigDefs.hpp \
  $(UTILS_H) \
  $(ANA_RTOP_BASE_H) \
  $(ANA_RTOP_BASE_C_H) \
  $(ANA_RTOP_OPS_LIB_H) \
  $(ANA_OPERATOR_VECTOR_INTERFACES_H)\
  $(ANA_OPERATOR_VECTOR_EXTENDED_INTERFACES_H)\
  $(ANA_OPERATOR_VECTOR_ANA_SUPPORT_H)\
	$(ANA_OPERATOR_VECTOR_ADAPTER_SUPPORT_H)\
  $(ANA_OPERATOR_VECTOR_ADAPTER_SUPPORT_SERIAL_H)\
  $(ANA_OPERATOR_VECTOR_ADAPTER_SUPPORT_MPI_H)

libthyra_a_SOURCES = \
  $(UTILS) \
  $(ANA_RTOP_BASE) \
  $(ANA_RTOP_BASE_C)


#We need the following to support separate build directories

# --------------------------------------------------------------------
# For using a special archiver
# --------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libthyra_a_AR = $(ALTERNATE_AR)

else

libthyra_a_AR = $(AR) cru

endif

# --------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# --------------------------------------------------------------------

# I have not been able to find any kind of "or" for Automake
# conditionals. JW

 libthyra_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# --------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# --------------------------------------------------------------------

lib_LIBRARIES = $(THYRA_LIB)

# --------------------------------------------------------------------
# List of all headers to install in $(includedir)
# --------------------------------------------------------------------

include_HEADERS = $(THYRA_H) 

# --------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# --------------------------------------------------------------------

# SPECIAL NOTE: ThyraConfig.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = Thyra_Config.h $(THYRA_H)

# --------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# --------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Library Target for Circular Dependencies
# ------------------------------------------------------------------------
all-libs:
	$(MAKE) libthyra.a
	cd $(top_builddir)/../teuchos ; make all-libs