// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//            PyTrilinos.Epetra: Python Interface to Epetra
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

%define DOCSTRING
"The Epetra module allows access to The Trilinos package Epetra.  Note
that the 'Epetra_' prefix has been stripped from all Epetra objects,
but that if imported with 'from PyTrilinos import Epetra', these
objects exist in the 'Epetra' python namespace.  Use the python help()
facility for local documentation on classes and methods, or see the
on-line documentation for more in-depth information."
%enddef

%module(package="PyTrilinos", docstring=DOCSTRING) Epetra

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

// Epetra includes
#include "Epetra_Version.h"
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_DistObject.h"
#include "Epetra_CompObject.h"
#include "Epetra_BLAS.h"
#include "Epetra_LAPACK.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_MapColoring.h"
#include "Epetra_DataAccess.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialDenseOperator.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialSymDenseMatrix.h"

// Local includes
#include "Epetra_NumPyVector.h"
#include "Epetra_VectorHelper.h"
#include "NumPyArray.h"
#include "NumPyWrapper.h"
%}

// Ignore directives
%ignore operator<<(ostream &, const Epetra_Object &);// From python, use __str__
%ignore Epetra_Object::Print(ostream &) const;
%ignore Epetra_SerialComm::operator=(const Epetra_SerialComm &);
%ignore Epetra_CompObject::operator=(const Epetra_CompObject &);
%ignore Epetra_CompObject::UpdateFlops(int) const;   // Use long int version
%ignore Epetra_CompObject::UpdateFlops(float) const; // Use double version
%ignore Epetra_BlockMap::operator=(const Epetra_BlockMap &);
%ignore Epetra_Map::operator=(const Epetra_Map &);
%ignore Epetra_LocalMap::operator=(const Epetra_LocalMap &);
%ignore Epetra_MultiVector::operator=(const Epetra_MultiVector &);
%ignore Epetra_MultiVector::operator[](int);         // See %extend MultiVector
%ignore Epetra_MultiVector::operator[](int) const;   //       __getitem__()
%ignore Epetra_MultiVector::operator()(int) const;
%ignore Epetra_IntVector::operator=(const Epetra_IntVector &);
%ignore Epetra_IntVector::operator[](int);           // See %extend IntVector
%ignore Epetra_IntVector::operator[](int) const;     //       __getitem__()
%ignore Epetra_CrsGraph::operator=(const Epetra_CrsGraph &);
%ignore Epetra_CrsGraph::operator[](int);            // See %extend CrsGraph
%ignore Epetra_CrsGraph::operator[](int) const;      //       __getitem__()
%ignore Epetra_MapColoring::operator[](int);         // See %extend Mapcoloring
%ignore Epetra_MapColoring::operator[](int) const;   //       __getitem__()
%ignore Epetra_MapColoring::operator()(int) const;
%ignore Epetra_CrsMatrix::operator=(const Epetra_CrsMatrix &);
%ignore Epetra_CrsMatrix::operator[](int);           // See %extend CrsMatrix
%ignore Epetra_CrsMatrix::operator[](int) const;     //       __getitem__()
// Epetra_VbrMatrix member function
// Solve(bool,bool,bool,Epetra_Vector const&,Epetra_Vector&) const doesn't
// appear to be implemented.  Apparently overridden by
// Solve(bool,bool,bool,Epetra_MultiVector const&,Epetra_MultiVector&) const
%ignore Epetra_VbrMatrix::operator=(const Epetra_VbrMatrix &);
%ignore Epetra_VbrMatrix::Solve(bool, bool, bool,
				Epetra_Vector const&, Epetra_Vector&) const;
%ignore Epetra_SerialDenseMatrix::operator=(const Epetra_SerialDenseMatrix &);
%ignore Epetra_SerialDenseMatrix::operator[](int);
%ignore Epetra_SerialDenseMatrix::operator[](int) const;
%ignore Epetra_SerialDenseMatrix::operator()(int,int) const;
%ignore Epetra_SerialDenseMatrix::A() const;
%ignore Epetra_SerialDenseVector::operator=(const Epetra_SerialDenseVector &);
%ignore Epetra_SerialDenseVector::operator[](int);
%ignore Epetra_SerialDenseVector::operator[](int) const;
%ignore Epetra_SerialDenseVector::operator()(int);
%ignore Epetra_SerialDenseVector::operator()(int) const;
%ignore NumPyArrayBase::print(std::ostream &) const;       // faciltated by __str__
%ignore NumPyArray::print(std::ostream &) const;           // faciltated by __str__
%ignore NumPyArrayContiguous::print(std::ostream &) const; // faciltated by __str__
%ignore NumPyArrayBase::getDataArray() const;
%ignore NumPyArrayBase::getArrayObject() const;

// Rename directives
%rename(Version             ) Epetra_Version;
%rename(Object              ) Epetra_Object;
%rename(Comm                ) Epetra_Comm;
%rename(SerialComm          ) Epetra_SerialComm;
%rename(BlockMap            ) Epetra_BlockMap;
%rename(Map                 ) Epetra_Map;
%rename(LocalMap            ) Epetra_LocalMap;
%rename(SrcDistObject       ) Epetra_SrcDistObject;
%rename(DistObject          ) Epetra_DistObject;
%rename(CompObject          ) Epetra_CompObject;
%rename(BLAS                ) Epetra_BLAS;
%rename(LAPACK              ) Epetra_LAPACK;
%rename(MultiVector         ) Epetra_MultiVector;
%rename(IntVector           ) Epetra_IntVector;
%rename(CrsGraph            ) Epetra_CrsGraph;
%rename(MapColoring         ) Epetra_MapColoring;
%rename(Operator            ) Epetra_Operator;
%rename(RowMatrix           ) Epetra_RowMatrix;
%rename(CrsMatrix           ) Epetra_CrsMatrix;
%rename(VbrMatrix           ) Epetra_VbrMatrix;
%rename(SerialDenseSolver   ) Epetra_SerialDenseSolver;
%rename(SerialDenseOperator ) Epetra_SerialDenseOperator;
%rename(SerialDenseMatrix   ) Epetra_SerialDenseMatrix;
%rename(SerialDenseVector   ) Epetra_SerialDenseVector;
%rename(SerialSymDenseMatrix) Epetra_SerialSymDenseMatrix;
%rename(NumPyVector         ) Epetra_NumPyVector;

// Auto-documentation feature
%feature("autodoc", "1");

// Typemap directives
%typemap(in) (int * Indices)
{
  PyArrayObject * array = NumPyWrapper::contiguous_typed_array($input, PyArray_INT,
							       0, 0);
  assert (array && "Function should have already checked this.");
  $1 = (int *) array->data;
}

%typemap(argout) (int & NumIndices, int *& Indices)
{
  // Decrement the current result object
  Py_XDECREF($result);
  // Check for bad result
  if (result == -1)
    PyErr_SetString(PyExc_RuntimeError, "Invalid row index");
  if (result == -2)
    PyErr_SetString(PyExc_RuntimeError, "Graph not completed");
  // Copy the indices into a python tuple
  $result = PyTuple_New(*$1);
  for (int i=0; i<*$1; ++i)
    PyTuple_SetItem($result, i, PyInt_FromLong((long) *$2[i])); 
}

// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"


// Epetra interface includes
using namespace std;
%include "Epetra_Version.h"
%include "Epetra_Object.h"
%include "Epetra_Comm.h"
%include "Epetra_SerialComm.h"
%include "Epetra_BlockMap.h"
%include "Epetra_Map.h"
%include "Epetra_LocalMap.h"
%include "Epetra_SrcDistObject.h"
%include "Epetra_DistObject.h"
%include "Epetra_CompObject.h"
%include "Epetra_BLAS.h"
%include "Epetra_LAPACK.h"
%include "Epetra_MultiVector.h"
%include "Epetra_Vector.h"
%include "Epetra_IntVector.h"
%include "Epetra_CrsGraph.h"
%include "Epetra_MapColoring.h"
%include "Epetra_Operator.h"
%include "Epetra_RowMatrix.h"
%include "Epetra_CrsMatrix.h"
%include "Epetra_VbrMatrix.h"
%include "Epetra_DataAccess.h"
%include "Epetra_SerialDenseSolver.h"
%include "Epetra_SerialDenseOperator.h"
%include "Epetra_SerialDenseMatrix.h"
%include "Epetra_SerialDenseVector.h"
%include "Epetra_SerialSymDenseMatrix.h"

// Local interface includes
%include "NumPyArray.h"
%include "Epetra_NumPyVector.h"

// Extensions
%extend Epetra_Object {
  string __str__() {
    stringstream os;
    self->Print(os);                  // Put the output in os
    string s = os.str();              // Extract the string from os
    return s.substr(0,s.length()-1);  // Return the string minus trailing \n
  }
}

%extend Epetra_MultiVector {
  double * & __getitem__(int i) {
    return self->operator[](i);
  }

  void Print() {
    self->Print(cout);
  }

  PyObject * Norm1() {
    int n = self->NumVectors();
    double result[n];
    int numVectors[1] = {n};
    int status        = self->Norm1(result);
    PyObject * output;
    if (n == 1) {
      output = Py_BuildValue("(id)", status, result[0]);
    } else {
      output = Py_BuildValue("(iO)", status, PyArray_FromDimsAndData(1,numVectors, PyArray_DOUBLE,
								     (char *)result));
    }
    return output;
  }

  PyObject * Norm2() {
    int n = self->NumVectors();
    double result[n];
    int numVectors[1] = {n};
    int status        = self->Norm2(result);
    PyObject * output;
    if (n == 1) {
      output = Py_BuildValue("(id)", status, result[0]);
    } else {
      output = Py_BuildValue("(iO)", status, PyArray_FromDimsAndData(1,numVectors, PyArray_DOUBLE,
								     (char *)result));
    }
    return output;
  }

  PyObject * NormInf() {
    int n = self->NumVectors();
    double result[n];
    int numVectors[1] = {n};
    int status        = self->NormInf(result);
    PyObject * output;
    if (n == 1) {
      output = Py_BuildValue("(id)", status, result[0]);
    } else {
      output = Py_BuildValue("(iO)", status, PyArray_FromDimsAndData(1,numVectors, PyArray_DOUBLE,
								     (char *)result));
    }
    return output;
  }
}

%extend Epetra_IntVector {
  int & __getitem__(int i) {
    return self->operator[](i);
  }
}

%extend Epetra_CrsGraph {
  int * __getitem__(int i) {
    return self->operator[](i);
  }
}

%extend Epetra_MapColoring {
  int & __getitem__(int i) {
    return self->operator[](i);
  }
}

%extend Epetra_CrsMatrix {
  double * __getitem__(int i) {
    return self->operator[](i);
  }
}

%extend Epetra_SerialDenseMatrix {

  double * __getitem__(int i) {
    return self->operator[](i);
  }

  PyObject * __getitem__(PyObject * args) {
    int i, j;
    if (!PyArg_ParseTuple(args, "ii", &i, &j)) {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return NULL;
    }
    double * column = self->operator[](j);
    return PyFloat_FromDouble(column[i]);
  }

  PyObject * __setitem__(PyObject * args, double val) {
    int i, j;
    if (!PyArg_ParseTuple(args, "ii", &i, &j)) {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return NULL;
    }
    double * column = self->operator[](j);
    column[i] = val;
    Py_INCREF(Py_None);
    return Py_None;
  }

}

%extend Epetra_SerialDenseVector {

  double __call__(int i) {
    return self->operator()(i);
  }

  double __getitem__(int i) {
    return self->operator[](i);
  }

  void __setitem__(int i, const double val) {
    double * column = self->Values();
    column[i] = val;
  }
}

%extend NumPyArrayBase {
  using namespace std;
  string __str__() {
    stringstream os;
    self->print(os);                  // Put the output in os
    string s = os.str();              // Extract the string from os
    return s.substr(0,s.length()-1);  // Return the string minus trailing \n
  }
}

// Python code.  Here we declare classes that inherit both from Epetra
// objects and UserArrays, to give users additional functionality.
%pythoncode %{

from UserArray import *

class Vector(UserArray,NumPyVector):
    def __init__(self, *args):
        NumPyVector.__init__(self, *args)
        UserArray.__init__(self,self.getArray(),'d',0,1)
    def __str__(self):
        return str(self.array)
%}
