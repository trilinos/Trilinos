# @HEADER
# ************************************************************************
# 
# RTOp: Interfaces and Support Software for Vector Reduction Transformation
#       Operations
#                Copyright (2006) Sandia Corporation
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
# Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package rtop subdirectory src built successfully."
	@echo ""

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

RTOP_INTERFACES_H = \
  $(srcdir)/interfaces/RTOpPack_Types.hpp\
  $(srcdir)/interfaces/RTOpPack_RTOpTDecl.hpp\
  $(srcdir)/interfaces/RTOpPack_RTOpT.hpp

RTOP_SUPPORT_H = \
  $(srcdir)/support/RTOpPack_MPI_apply_op.hpp\
  $(srcdir)/support/RTOpPack_RTOpServer.hpp\
  $(srcdir)/support/RTOpPack_RTOpServerDecl.hpp\
  $(srcdir)/support/RTOpPack_RTOpTHelpers.hpp\
  $(srcdir)/support/RTOpPack_SPMD_apply_op_decl.hpp\
  $(srcdir)/support/RTOpPack_TOpUnaryFuncPtrDecl.hpp\
  $(srcdir)/support/RTOpPack_MPI_apply_op_decl.hpp\
  $(srcdir)/support/RTOpPack_TOpUnaryFuncPtr.hpp

RTOP_SUPPORT = \
  $(srcdir)/support/RTOpPack_MPI_apply_op.cpp\
  $(srcdir)/support/RTOpPack_SPMD_apply_op.cpp

RTOP_SUPPORT_C_H = \
  $(srcdir)/support/RTOp_MPI_config.h\
  $(srcdir)/support/RTOp_mpi.h\
  $(srcdir)/support/RTOp_parallel_helpers.h

RTOP_SUPPORT_C = \
  $(srcdir)/support/RTOp_mpi.c\
  $(srcdir)/support/RTOp_parallel_helpers.c

RTOP_OPS_LIB_H = \
  $(srcdir)/ops_lib/RTOpPack_ROpCountNanInf.hpp\
  $(srcdir)/ops_lib/RTOpPack_ROpNorm1.hpp\
  $(srcdir)/ops_lib/RTOpPack_TOpAssignVectors.hpp\
  $(srcdir)/ops_lib/RTOpPack_ROpDotProd.hpp\
  $(srcdir)/ops_lib/RTOpPack_ROpNorm2.hpp\
  $(srcdir)/ops_lib/RTOpPack_TOpEleWiseDivide.hpp\
  $(srcdir)/ops_lib/RTOpPack_ROpGetSubVector.hpp\
  $(srcdir)/ops_lib/RTOpPack_ROpNormInf.hpp\
  $(srcdir)/ops_lib/RTOpPack_TOpEleWiseProd.hpp\
  $(srcdir)/ops_lib/RTOpPack_TOpEleWiseProdUpdate.hpp\
  $(srcdir)/ops_lib/RTOpPack_ROpMax.hpp\
  $(srcdir)/ops_lib/RTOpPack_ROpSum.hpp\
  $(srcdir)/ops_lib/RTOpPack_TOpLinearCombination.hpp\
  $(srcdir)/ops_lib/RTOpPack_ROpMaxIndex.hpp\
  $(srcdir)/ops_lib/RTOpPack_ROpWeightedNorm2.hpp\
  $(srcdir)/ops_lib/RTOpPack_TOpRandomize.hpp\
  $(srcdir)/ops_lib/RTOpPack_ROpMaxIndexLessThanBound.hpp\
  $(srcdir)/ops_lib/RTOpPack_TOpAXPY.hpp\
  $(srcdir)/ops_lib/RTOpPack_TOpReciprocal.hpp\
  $(srcdir)/ops_lib/RTOpPack_ROpMin.hpp\
  $(srcdir)/ops_lib/RTOpPack_TOpAbs.hpp\
  $(srcdir)/ops_lib/RTOpPack_TOpScaleVector.hpp\
  $(srcdir)/ops_lib/RTOpPack_ROpMinIndex.hpp\
  $(srcdir)/ops_lib/RTOpPack_TOpAddScalar.hpp\
  $(srcdir)/ops_lib/RTOpPack_TOpSetSubVector.hpp\
  $(srcdir)/ops_lib/RTOpPack_ROpMinIndexGreaterThanBound.hpp\
  $(srcdir)/ops_lib/RTOpPack_TOpAssignScalar.hpp\
  $(srcdir)/ops_lib/RTOpPack_SUNDIALS_Ops.hpp

#Specify locations of include files within Trilinos

include $(top_builddir)/Makefile.export.rtop

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(RTOP_INCLUDES))
else
EXPORT_INCLUDES = $(RTOP_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# --------------------------------------------------------------------
# RTOp library specifications
# --------------------------------------------------------------------

RTOP_LIB = librtop.a

RTOP_H = \
  $(srcdir)/RTOp_ConfigDefs.hpp \
  $(RTOP_INTERFACES_H) \
  $(RTOP_SUPPORT_H) \
  $(RTOP_SUPPORT_C_H) \
  $(RTOP_OPS_LIB_H)

librtop_a_SOURCES = \
  $(RTOP_SUPPORT_C) \
  $(RTOP_SUPPORT)

#We need the following to support separate build directories

# --------------------------------------------------------------------
# For using a special archiver
# --------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

librtop_a_AR = $(ALTERNATE_AR)

else

librtop_a_AR = $(AR) cru

endif

# --------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# --------------------------------------------------------------------

# I have not been able to find any kind of "or" for Automake
# conditionals. JW

 librtop_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# --------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# --------------------------------------------------------------------

lib_LIBRARIES = $(RTOP_LIB)

# --------------------------------------------------------------------
# List of all headers to install in $(includedir)
# --------------------------------------------------------------------

include_HEADERS = $(RTOP_H) 

# --------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# --------------------------------------------------------------------

# SPECIAL NOTE: RTOpConfig.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = RTOp_Config.h $(RTOP_H)

# --------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# --------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Library Target for Circular Dependencies
# ------------------------------------------------------------------------
all-libs:
	$(MAKE) librtop.a
	cd $(top_builddir)/../teuchos ; $(MAKE) all-libs
