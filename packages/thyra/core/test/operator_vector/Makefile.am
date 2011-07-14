#@HEADER
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
#@HEADER

EXEEXT = .exe

#
# Actual programs
#

TEST_PRODUCT_SPACE = test_product_space
TEST_SCALAR_PRODUCT = test_scalar_product
TEST_COMPOSITE_LINEAR_OPS = test_composite_linear_ops
TEST_STD_OPS = test_std_ops
TEST_HANDLES = test_handles
TEST_BLOCK_OP = test_block_op
TEST_LINEAR_COMBINATION = test_linear_combination

noinst_PROGRAMS = $(TEST_PRODUCT_SPACE) $(TEST_SCALAR_PRODUCT) \
$(TEST_COMPOSITE_LINEAR_OPS) $(TEST_STD_OPS) $(TEST_HANDLES) \
$(TEST_BLOCK_OP) $(TEST_LINEAR_COMBINATION)

test_product_space_SOURCES = $(srcdir)/test_product_space.cpp
test_product_space_DEPENDENCIES = $(common_dependencies)
test_product_space_LDADD = $(common_ldadd)

test_scalar_product_SOURCES = $(srcdir)/test_scalar_product.cpp
test_scalar_product_DEPENDENCIES = $(common_dependencies)
test_scalar_product_LDADD = $(common_ldadd)

test_composite_linear_ops_SOURCES = $(srcdir)/test_composite_linear_ops.cpp
test_composite_linear_ops_DEPENDENCIES = $(common_dependencies)
test_composite_linear_ops_LDADD = $(common_ldadd)

test_std_ops_SOURCES = $(srcdir)/test_std_ops.cpp
test_std_ops_DEPENDENCIES = $(common_dependencies)
test_std_ops_LDADD = $(common_ldadd)

test_handles_SOURCES = $(srcdir)/test_handles.cpp
test_handles_DEPENDENCIES = $(common_dependencies)
test_handles_LDADD = $(common_ldadd)

test_block_op_SOURCES = $(srcdir)/test_block_op.cpp
test_block_op_DEPENDENCIES = $(common_dependencies)
test_block_op_LDADD = $(common_ldadd)

test_linear_combination_SOURCES = $(srcdir)/test_linear_combination.cpp
test_linear_combination_DEPENDENCIES = $(common_dependencies)
test_linear_combination_LDADD = $(common_ldadd)

#
# Common stuff
#

common_dependencies = \
  $(top_builddir)/src/libthyra.a

include $(top_builddir)/Makefile.export.thyra

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(THYRA_INCLUDES))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(THYRA_LIBS))
else
EXPORT_INCLUDES = $(THYRA_INCLUDES)
EXPORT_LIBS = $(THYRA_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

common_ldadd = $(EXPORT_LIBS)

# This is already added part of THYRA_LIBS and therefore automake does not need to add this again!
LIBS =

# Other files to add to the tarball
EXTRA_DIST = compare-thyra-spmd-multivector-files.py
