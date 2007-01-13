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

%define EPETRAEXT_DOCSTRING
"The EpetraExt module allows access to The Trilinos package EpetraExt.
Use the python help() facility for local documentation on classes and
methods, or see the on-line documentation for more in-depth
information.

The most important classes of the EpetraExt module are:
*) Graph coloring classes:
   - MapColoring
   - MapColoringIndex
*) Input functions:
   - BlockMapToMatrixMarketFile()
   - MatrixMarketFileToMultiVector()
   - MatrixMarketFileToBlockMap()
   - MatrixMarketFileToMap()
   - MatrixMarketFileToCrsMatrix()
*) Output functions:
   - BlockMapToMatrixMarketFile()
   - RowMatrixToMatrixMarketFile()
   - MultiVectorToMatrixMarketFile()
*) Input/Output classes:
   - HDF5
   - XML
*) Matrix-Matrix functions:
   - Add()
   - Multiply()
"
%enddef

%module(package   = "PyTrilinos",
	directors = "1",
	docstring = EPETRAEXT_DOCSTRING) EpetraExt

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
#include "EpetraExt_XMLWriter.h" // FIXME: memory management still scary...

%}

// Ignore directives
%ignore EpetraExt::HDF5::Read;
%ignore EpetraExt::XMLReader::Read;

using namespace std;
// C++ STL support
%include "stl.i"

// Epetra interface import
%import "Teuchos.i"
%import "Epetra.i"

// Typemaps
%define OUTPUT_ARGUMENT(ClassName)
%typemap(in,numinputs=0) ClassName *& (ClassName * _object) {
  $1 = &_object;
}
%typemap(argout) ClassName *& {
  PyObject * obj1;
  PyObject * obj2;
  static swig_type_info * swig_CN_ptr = SWIG_TypeQuery("ClassName *");
  obj1 = SWIG_NewPointerObj((void*)(*$1), swig_CN_ptr, 1);
  if (result < 0) obj1 = Py_BuildValue("");
  if (!PyTuple_Check($result)) $result = Py_BuildValue("(O)", $result);
  obj2 = Py_BuildValue("(O)", obj1);
  $result = PySequence_Concat($result,obj2);
}
%enddef

%define OUTPUT_EPETRA_ARRAY_ARGUMENT(ClassName)
%typemap(in,numinputs=0) Epetra_ ## ClassName *& (Epetra_ ## ClassName * _object) {
  $1 = &_object;
}
%typemap(argout) Epetra_ ## ClassName *& {
  PyObject * obj1;
  PyObject * obj2;
  static swig_type_info * swig_NP_ptr = SWIG_TypeQuery("Epetra_NumPy" "ClassName *");
  Epetra_NumPy ## ClassName * npa = new Epetra_NumPy ## ClassName(**$1);
  obj1 = SWIG_NewPointerObj((void*)npa, swig_NP_ptr, 1);
  if (result < 0) obj1 = Py_BuildValue("");
  if (!PyTuple_Check($result)) $result = Py_BuildValue("(O)", $result);
  obj2 = Py_BuildValue("(O)", obj1);
  $result = PySequence_Concat($result,obj2);
}
%enddef

OUTPUT_ARGUMENT(Epetra_BlockMap )
OUTPUT_ARGUMENT(Epetra_Map      )
OUTPUT_ARGUMENT(Epetra_CrsMatrix)

OUTPUT_EPETRA_ARRAY_ARGUMENT(MultiVector)
OUTPUT_EPETRA_ARRAY_ARGUMENT(Vector     )

// EpetraExt interface includes
%include "EpetraExt_Version.h"
%include "EpetraExt_HDF5.h"
%include "EpetraExt_XMLReader.h"
%include "EpetraExt_XMLWriter.h"

%include "EpetraExt_Transform.h"
%template () std::vector<Epetra_IntVector>;
%template () EpetraExt::Transform<Epetra_CrsGraph, Epetra_MapColoring>;
%template () EpetraExt::Transform<Epetra_CrsGraph, std::vector<Epetra_IntVector,
							       std::allocator<Epetra_IntVector> > >;
%template () EpetraExt::StructuralTransform<Epetra_CrsGraph, Epetra_MapColoring>;
%template () EpetraExt::StructuralTransform<Epetra_CrsGraph, std::vector<Epetra_IntVector> >;

%include "EpetraExt_MapColoring.h"
%include "EpetraExt_MapColoringIndex.h"

%include "EpetraExt_MultiVectorIn.h"
%include "EpetraExt_MultiVectorOut.h"
%include "EpetraExt_CrsMatrixIn.h"
%include "EpetraExt_RowMatrixOut.h"
%include "EpetraExt_BlockMapIn.h"
%include "EpetraExt_BlockMapOut.h"

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
%define READ_EPETRA_CLASS(ClassName)
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

namespace EpetraExt {

#ifdef HAVE_EPETRAEXT_HDF5
  %extend HDF5 {

    READ_EPETRA_CLASS(BlockMap   )
    READ_EPETRA_CLASS(Map        )
    READ_EPETRA_CLASS(MultiVector)
    READ_EPETRA_CLASS(CrsGraph   )
    READ_EPETRA_CLASS(CrsMatrix  )
    READ_EPETRA_CLASS(IntVector  )
  }    // HDF5
#endif

  %extend XMLReader {

    READ_EPETRA_CLASS(Map        )
    READ_EPETRA_CLASS(MultiVector)
    READ_EPETRA_CLASS(CrsGraph   )
    READ_EPETRA_CLASS(CrsMatrix  )
  }    // XMLReader

}    // EpetraExt


%inline %{
  namespace EpetraExt
  {
    int Add(Epetra_CrsMatrix& A, const bool flag, const double ValA,
            Epetra_CrsMatrix& B, const double ValB)
    {
      EpetraExt::MatrixMatrix M;
      return(M.Add(A, flag, ValA, B, ValB));
    }
  }
%}

// Python code.
%pythoncode %{
  __version__ = EpetraExt_Version().split()[2]
%}
