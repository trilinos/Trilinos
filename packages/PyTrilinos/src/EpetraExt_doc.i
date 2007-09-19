// -*- C++ -*-
// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%feature("autodoc",
"MatlabFileToCrsMatrix(str filename, Epetra.Comm) -> Epetra.CrsMatrix

Return a CrsMatrix read from a matlab file.")
EpetraExt::MatlabFileToCrsMatrix;

%feature("autodoc",
"MatrixMarketFileToBlockMap(str filename, Epetra.Comm) -> Epetra.BlockMap

Return a BlockMap read from a matrix market file.")
EpetraExt::MatrixMarketFileToBlockMap;

%feature("autodoc",
"MatrixMarketFileToBlockMaps(str filename, Epetra.Comm) ->
    (Epetra.BlockMap rowMap, Epetra.BlockMap colMap, Epetra.BlockMap rangeMap,
     Epetra.BlockMap domainMap)

Return a tuple of BlockMaps read from a matrix market file.  The
BlockMaps, listed in order, are the row map, the column map, the range
map and the domain map.")
EpetraExt::MatrixMarketFileToBlockMaps;

%feature("autodoc",
"MatrixMarketFileToCrsMatrix(str filename, Epetra.Map rowMap, Epetra.Map
    colMap=None, Epetra.Map rangeMap=None, Epetra.Map domainMap=None) ->
    Epetra.CrsMatrix

Return a CrsMatrix read from a matrix market file.")
EpetraExt::MatrixMarketFileToCrsMatrix;

%feature("autodoc",
"MatrixMarketFileToMap(str filename, Epetra.Comm) -> Epetra.Map

Return a Map read from a matrix market file.")
EpetraExt::MatrixMarketFileToMap;

%feature("autodoc",
"MatrixMarketFileToMultiVector(str filename, Epetra.BlockMap) ->
    Epetra.MultiVector

Return a MultiVector read from a matix market file.")
EpetraExt::MatrixMarketFileToMultiVector;

// EpetraExt.HDF5 class

%feature("docstring")
EpetraExt::HDF5::ReadBlockMap
"
Return a BlockMap read from an HDF5 file specified by filename 'name'.
"

%feature("docstring")
EpetraExt::HDF5::ReadMap
"
Return a Map read from an HDF5 file specified by filename 'name'.
"

%feature("docstring")
EpetraExt::HDF5::ReadIntVector
"
Return an IntVector read from an HDF5 file specified by filename
'name'.
"

%feature("docstring")
EpetraExt::HDF5::ReadMultiVector
"
Return a MultiVector read from an HDF5 file specified by filename
'name'.
"

%feature("docstring")
EpetraExt::HDF5::ReadCrsGraph
"
Return a CrsGraph read from an HDF5 file specified by filename 'name'.
"

%feature("docstring")
EpetraExt::HDF5::ReadCrsMatrix
"
Return a CrsMatrix read from an HDF5 file specified by filename
'name'.
"

// EpetraExt::XMLReader class

%feature("docstring")
EpetraExt::XMLReader::ReadMap
"
Return a Map read from an XML file specified by filename 'name'.
"

%feature("docstring")
EpetraExt::XMLReader::ReadMultiVector
"
Return a MultiVector read from an XML file specified by filename
'name'.
"

%feature("docstring")
EpetraExt::XMLReader::ReadCrsGraph
"
Return a CrsGraph read from an XML file specified by filename 'name'.
"

%feature("docstring")
EpetraExt::XMLReader::ReadCrsMatrix
"
Return a CrsMatrix read from an XML file specified by filename 'name'.
"

// EpetraExt::Add function

%feature("autodoc",
"Add(Epetra.CrsMatrix A, bool flag, float valA, Epetra.CrsMatrix B,
    float valB) -> int

Compute B <- valA * A + valB * B.  If flag is True, use the transpose
of A.  B must either have the structure of A+B or not yet have
FillComplete() called on it.")
EpetraExt::Add;

// EpetraExt::Multiply function

%feature("autodoc",
"Multiply(Epetra.CrsMatrix A, bool transposeA, Epetra.CrsMatrix B, bool
    transposeB, Epetra.CrsMatrix C) -> int

Compute C <- A * B, where transposeA and transposeB control the
transposition of A and B respectively.  C must have the structure of A
* B, or not yet have FillComplete() called on it.")
EpetraExt::Multiply;
