#@HEADER
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
#@HEADER

# The following line helps the test harness recover from build errors.

all-local:

SUBDIRS = IlukGraph               \
          SetParameters           \
          CompareWithAztecOO      \
          Container               \
          OverlappingRowMatrix    \
          PreconditionerFactory   \
          Relaxation              \
          ILU                     \
          IC                      \
	  UseTranspose            \
	  BlockCheby              \
	  TestAll

# CrsRiluk 

if HAVE_IFPACK_AZTECOO
IFPACK_COMPARE_WITH_AZTECOO_TEST = CompareWithAztecOO/TestCompareWithAztecOO.exe
endif

if HAVE_IFPACK_AMESOS
IFPACK_CONTAINER_TEST = Container/TestContainer.exe
endif

if HAVE_IFPACK_AZTECOO
if HAVE_IFPACK_AMESOS
IFPACK_OVERLAPPING_ROW_MATRIX_TEST = OverlappingRowMatrix/TestOverlappingRowMatrix.exe
endif
endif

if HAVE_IFPACK_AZTECOO
IFPACK_PRECONDITIONER_FACTORY_TEST = PreconditionerFactory/TestPreconditionerFactory.exe
endif

if HAVE_IFPACK_AZTECOO
IFPACK_RELAXATION_TEST = Relaxation/TestRelaxation.exe BlockCheby/TestBlockCheby.exe
endif

if HAVE_IFPACK_AMESOS
IFPACK_ILU_TEST = ILU/TestILU.exe
endif

if HAVE_IFPACK_AMESOS
IFPACK_IC_TEST = IC/TestIC.exe
endif

TESTS = \
CrsRiluk/CrsRiluk_test.exe \
IlukGraph/IlukGraph_test.exe \
SetParameters/SetParameters.exe \
UseTranspose/TestUseTranspose.exe \
$(IFPACK_COMPARE_WITH_AZTECOO_TEST) \
$(IFPACK_CONTAINER_TEST) \
$(IFPACK_OVERLAPPING_ROW_MATRIX_TEST) \
$(IFPACK_PRECONDITIONER_FACTORY_TEST) \
$(IFPACK_RELAXATION_TEST) \
$(IFPACK_ILU_TEST) \
$(IFPACK_IC_TEST)

EXTRA_DIST = definition

