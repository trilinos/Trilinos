#@HEADER
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
#@HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package epetra subdirectory test built successfully."
	@echo ""

SUBDIRS = src BlockMap CrsGraph CrsMatrix CrsRectMatrix Directory \
FECrsGraph FECrsMatrix FEVbrMatrix FEVector ImportExport \
IntSerialDense Map MapColoring \
MultiVector Object RowMatrixTransposer SerialDense SerialSpdDense \
VbrMatrix Vector Comm BasicPerfTest

TESTS = \
BlockMap/BlockMap_test.exe \
Comm/Comm_test.exe \
CrsGraph/CrsGraph_test.exe \
CrsMatrix/CrsMatrix_test.exe \
CrsRectMatrix/CrsRectMatrix_test.exe \
Directory/Directory_test.exe \
FECrsGraph/FECrsGraph_test.exe \
FECrsMatrix/FECrsMatrix_test.exe \
FEVbrMatrix/FEVbrMatrix_test.exe \
FEVector/FEVector_test.exe \
ImportExport/ImportExport_test.exe \
IntSerialDense/IntSerialDense_test.exe \
Map/Map_test.exe \
MapColoring/MapColoring_test.exe \
MultiVector/MultiVector_test.exe \
Object/Object_test.exe \
RowMatrixTransposer/RowMatrixTransposer_test.exe \
SerialDense/SerialDense_test.exe \
SerialSpdDense/SerialSpdDense_test.exe \
VbrMatrix/VbrMatrix_test.exe \
Vector/Vector_test.exe

EXTRA_DIST = definition

