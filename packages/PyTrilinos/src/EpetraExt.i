// -*- c++ -*-

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

%define %epetraext_docstring
"
PyTrilinos.EpetraExt is the python interface to Trilinos package
EpetraExt:

    http://software.sandia.gov/trilinos/packages/epetraext

The purpose of EpetraExt is to provide various extensions to Epetra
that were not considered appropriate for the Epetra package.  These
extensions include I/O, matrix-matrix operations and graph coloring.

Currently, only a subset of EpetraExt classes and functions have
python interfaces, including the following user-level classes:

    * XMLReader                      - Read Epetra data from an XML file
    * XMLWriter                      - Write Epetra data as an XML file
    * CrsGraph_MapColoring           - Compute a graph coloring
    * CrsGraph_MapColoringIndex      - Compute indexes for a graph coloring

and functions:

    * MatrixMarketFileToBlockMap     - Read a BlockMap from an MM file
    * MatrixMarketFileToBlockMaps    - Read BlockMaps from an MM file
    * MatrixMarketFileToMap          - Read a Map from an MM file
    * MatrixMarketFileToMultiVector  - Read a MultiVector from an MM file
    * MatrixMarketFileToCrsMatrix    - Read a CrsMatrix from an MM file
    * MatlabFileToCrsMatrix          - Read a CrsMatrix from an ML file
    * BlockMapToMatrixMarketFile     - Write a BlockMap to an MM file
    * MultiVectorToMatrixMarketFile  - Write a MultiVector to an MM file
    * MultiVectorToMatlabFile        - Write a MultiVector to an ML file
    * RowMatrixToMatrixMarketFile    - Write a RowMatrix to an MM file
    * RowMatrixToMatlabFile          - Write a RowMatrix to an ML file
    * Add                            - Add two CrsMatrix objects

For examples of usage, please consult the following scripts in the
example subdirectory of the PyTrilinos package:

    * exEpetraExt_IO_MatrixMarket.py
    * exEpetraExt_IO_XML.py
    * exEpetraExt_MatrixMatrix.py
"
%enddef

%module(package   = "PyTrilinos",
	directors = "1",
	docstring = %epetraext_docstring) EpetraExt

%{
// System includes
#include <vector>

// Configuration includes
#include "PyTrilinos_config.h"

// Teuchos includes
#include "Teuchos_PythonParameter.h"

// Epetra includes
#include "Epetra_Object.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_DistObject.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_LocalMap.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MapColoring.h"

// Epetra python includes
#include "NumPyImporter.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyIntVector.h"

// EpetraExt includes
#include "EpetraExt_Version.h"
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_HDF5.h"
#include "EpetraExt_XMLReader.h"
#include "EpetraExt_XMLWriter.h"
%}

////////////
// Macros //
////////////

// The overloaded HDF5 and XMLReader Read() methods cannot be
// type-disambiguated in python.  We therefore replace selected
// overloaded Read() methods with a python version that has the type
// in the method name.  For example,
//
//   void HDF5::Read(std::string, Epetra_Map *&)
//
// is translated from C++ to python as
//
//   HDF5.ReadMap(str) -> Epetra.Map
//
// These translations are made possible by the following macro:
//
%define %epetraext_read_method(ClassName)
Epetra_ ## ClassName * Read ## ClassName(string name)
{
  Epetra_ ## ClassName * obj = NULL;
  try {
    self->Read(name, obj);
  } catch(std::logic_error e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    return NULL;
  }
  return obj;
}
%enddef

// C++ STL support
using namespace std;
%include "stl.i"

// Trilinos interface support
%import "Teuchos.i"
%import "Epetra.i"

///////////////////////////////
// EpetraExt_Version support //
///////////////////////////////
%include "EpetraExt_Version.h"
%pythoncode %{
__version__ = EpetraExt_Version().split()[2]
%}

////////////////////////////
// EpetraExt_HDF5 support //
////////////////////////////
%ignore EpetraExt::HDF5::Read;
%include "EpetraExt_HDF5.h"
namespace EpetraExt {
#ifdef HAVE_EPETRAEXT_HDF5
  %extend HDF5 {
    %epetraext_read_method(BlockMap   )
    %epetraext_read_method(Map        )
    %epetraext_read_method(MultiVector)
    %epetraext_read_method(CrsGraph   )
    %epetraext_read_method(CrsMatrix  )
    %epetraext_read_method(IntVector  )
  }    // HDF5
#endif
}

/////////////////////////////////
// EpetraExt_XMLReader support //
/////////////////////////////////
%ignore EpetraExt::XMLReader::Read;
%include "EpetraExt_XMLReader.h"
namespace EpetraExt {
  %extend XMLReader {
    %epetraext_read_method(Map        )
    %epetraext_read_method(MultiVector)
    %epetraext_read_method(CrsGraph   )
    %epetraext_read_method(CrsMatrix  )
  }    // XMLReader
}

/////////////////////////////////
// EpetraExt_XMLWriter support //
/////////////////////////////////
%include "EpetraExt_XMLWriter.h"

/////////////////////////////////
// EpetraExt_Transform support //
/////////////////////////////////
%include "EpetraExt_Transform.h"
%template () std::vector<Epetra_IntVector>;
%template () EpetraExt::Transform<Epetra_CrsGraph, Epetra_MapColoring>;
%template () EpetraExt::Transform<Epetra_CrsGraph, std::vector<Epetra_IntVector,
							       std::allocator<Epetra_IntVector> > >;
%template () EpetraExt::StructuralTransform<Epetra_CrsGraph, Epetra_MapColoring>;
%template () EpetraExt::StructuralTransform<Epetra_CrsGraph, std::vector<Epetra_IntVector> >;

///////////////////////////////
// EpetraExt_Version support //
///////////////////////////////
%include "EpetraExt_MapColoring.h"

///////////////////////////////
// EpetraExt_Version support //
///////////////////////////////
%include "EpetraExt_MapColoringIndex.h"

///////////////////////////////
// EpetraExt_Version support //
///////////////////////////////
%include "EpetraExt_MultiVectorIn.h"

///////////////////////////////
// EpetraExt_Version support //
///////////////////////////////
%include "EpetraExt_MultiVectorOut.h"

///////////////////////////////
// EpetraExt_Version support //
///////////////////////////////
%include "EpetraExt_CrsMatrixIn.h"

///////////////////////////////
// EpetraExt_Version support //
///////////////////////////////
%include "EpetraExt_RowMatrixOut.h"

///////////////////////////////
// EpetraExt_Version support //
///////////////////////////////
%include "EpetraExt_BlockMapIn.h"

///////////////////////////////
// EpetraExt_Version support //
///////////////////////////////
%include "EpetraExt_BlockMapOut.h"

/////////////////////////////
// EpetraExt.Add() support //
/////////////////////////////
%inline %{
  namespace EpetraExt {
    int Add(Epetra_CrsMatrix& A, const bool flag, const double ValA,
            Epetra_CrsMatrix& B, const double ValB) {
      EpetraExt::MatrixMatrix M;
      return(M.Add(A, flag, ValA, B, ValB));
    }
  }
%}
