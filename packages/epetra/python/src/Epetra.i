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

%define EPETRA_DOCSTRING
"The Epetra module allows access to The Trilinos package Epetra.  Note
that the 'Epetra_' prefix has been stripped from all Epetra objects,
but that if imported with 'from PyTrilinos import Epetra', these
objects exist in the 'Epetra' python namespace.  Use the python help()
facility for local documentation on classes and methods, or see the
on-line documentation for more in-depth information.


Brief Description
=================

The Epetra module offers a suite of serial and distributed linear algebra
objects, like vectors, multivectors, sparse matrices. Serial dense vectors and
matrices are also supported. Using Epetra, several BLAS and LAPACK routines
can be easily accessed through Python. The Epetra communicator is the basic
intra-processor class of Trilinos.

The most important classes of the Epetra module is:
- PyComm: basic serial or parallel communicator
- SerialDenseVector: serial vector, to be used with a serial dense matrix
- SerialDenseMatrix: serial matrix, allocated as dense
- SerialDenseProblem: to solve serial dense problems
- Map: layout of distributed objects across processors
- Vector: distributed vector
- MultiVector: series of Vector's with the same Map
- CrsMatrix: distributed matrix with row access
- LinearProblem: to solve distributed linear systems
- Time: utility class to track CPU times
- Import: to easily redistribute distributed objects


Example of usage: Creating a Distributed Matrix (1)
===================================================

The following example builds a distributed tridiagonal matrix. The
matrix has size 10 and the elements are distributed linearly among the
processors.  It can be run is serial or parallel, depending on how
Trilinos was configured.

from PyTrilinos import Epetra
Comm  = Epetra.PyComm()
NumGlobalRows = 10
Map   = Epetra.Map(NumGlobalRows, 0, Comm)
A     = Epetra.CrsMatrix(Epetra.Copy, Map, 0);
MyGlobalElements = Map.MyGlobalElements()
# Loop over all local rows to create a tridiagonal matrix
for i in MyGlobalElements:
  if i != 0:
    A[i, i - 1] = -1.0
  if i != NumGlobalRows - 1:
    A[i, i + 1] = -1.0
  A[i, i] = 2.0

# transform the matrix into local representation -- no entries
# can be added after this point. However, existing entries can be
# modified.
ierr = A.FillComplete();

# Now accessing the matrix nonzeros on each row
for i in MyGlobalElements:
  Indices, Values = A[i]
  print Indices, Values


Example of usage: Creating a Distributed Matrix (2)
===================================================

As in the previous example, but more efficiently.

from PyTrilinos import Epetra
Comm  = Epetra.PyComm()
NumGlobalRows = 10
Map   = Epetra.Map(NumGlobalRows, 0, Comm)
A     = Epetra.CrsMatrix(Epetra.Copy, Map, 0);
NumMyRows = Map.NumMyElements()
# Loop over all local rows to create a tridiagonal matrix
for ii in xrange(0, NumMyRows):
  # `i' is the global ID of local ID `ii'
  i = Map.GID(ii)
  if i != 0:
    Indices = [i, i - 1]
    Values = [2.0, -1.0]
  if i != NumGlobalRows - 1:
    Indices = [i, i + 1]
    Values = [2.0, -1.0]
  else:
    Indices = [i, i - 1, i + 1]
    Values = [2.0, -1.0, -1.0];
  A.InsertGlobalValues(i, Values, Indices);

# transform the matrix into local representation -- no entries
# can be added after this point. However, existing entries can be
# modified.
ierr = A.FillComplete();


Example of usage: Matrix-Vector product
=======================================

This example show how to multiply a matrix (diagonal here for simplicity) by a
vector.

from PyTrilinos import Epetra
Comm  = Epetra.PyComm()
NumGlobalRows = 10
Map   = Epetra.Map(NumGlobalRows, 0, Comm)
A     = Epetra.CrsMatrix(Epetra.Copy, Map, 0);
MyGlobalElements = Map.MyGlobalElements()
for i in MyGlobalElements:
  A[i, i] = 1.0 + i

A.FillComplete()

X = Epetra.Vector(Map); X.Random()
Y = Epetra.Vector(Map)

UseTranspose = False
A.Multiply(UseTranspose, X, Y)
  

Notes
=====

An Epetra.Vector object is at the same time a Numeric vector, and it
can therefore be used everywhere Numeric vectors are accepted.
"
%enddef

%module(package="PyTrilinos", directors="1", docstring=EPETRA_DOCSTRING) Epetra

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

#include "Epetra_ConfigDefs.h"

// Epetra includes
#include "Epetra_Version.h"
#include "Epetra_Object.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_DistObject.h"
#include "Epetra_CompObject.h"
#include "Epetra_BLAS.h"
#include "Epetra_LAPACK.h"
#include "Epetra_IntVector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_MapColoring.h"
#include "Epetra_DataAccess.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialDenseOperator.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialSymDenseMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_CombineMode.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Time.h"

// Local includes
#include "FileStream.h"
#include "NumPyArray.h"
#include "NumPyWrapper.h"
#include "PyEpetra_Utils.h"  

#include "Epetra_PyOperator.h"
#include "Epetra_PyRowMatrix.h"
%}

%feature("director") PyOperator;
%feature("director") PyRowMatrix;

// Ignore directives
%ignore operator<<(ostream &, const Epetra_Object &);// From python, use __str__
%ignore Epetra_Object::Print(ostream &) const;
%ignore Epetra_CompObject::operator=(const Epetra_CompObject &);
%ignore Epetra_CompObject::UpdateFlops(int) const;   // Use long int version
%ignore Epetra_CompObject::UpdateFlops(float) const; // Use double version
%ignore Epetra_BlockMap::operator=(const Epetra_BlockMap &);
%ignore Epetra_Map::operator=(const Epetra_Map &);
%ignore Epetra_LocalMap::operator=(const Epetra_LocalMap &);
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
%ignore Epetra_FECrsMatrix::operator=(const Epetra_FECrsMatrix &);
%ignore Epetra_FECrsMatrix::operator[](int);           // See %extend FECrsMatrix
%ignore Epetra_FECrsMatrix::operator[](int) const;     //       __getitem__()
// Epetra_VbrMatrix member function
// Solve(bool,bool,bool,Epetra_Vector const&,Epetra_Vector&) const does not
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
%ignore Epetra_IntSerialDenseMatrix::operator=(const Epetra_IntSerialDenseMatrix &);
%ignore Epetra_IntSerialDenseMatrix::operator[](int);
%ignore Epetra_IntSerialDenseMatrix::operator[](int) const;
%ignore Epetra_IntSerialDenseMatrix::operator()(int,int) const;
%ignore Epetra_IntSerialDenseMatrix::A() const;
%ignore Epetra_IntSerialDenseVector::operator=(const Epetra_IntSerialDenseVector &);
%ignore Epetra_IntSerialDenseVector::operator[](int);
%ignore Epetra_IntSerialDenseVector::operator[](int) const;
%ignore Epetra_IntSerialDenseVector::operator()(int);
%ignore Epetra_IntSerialDenseVector::operator()(int) const;
%ignore NumPyArrayBase::print(std::ostream &) const;       // faciltated by __str__
%ignore NumPyArray::print(std::ostream &) const;           // faciltated by __str__
%ignore NumPyArrayContiguous::print(std::ostream &) const; // faciltated by __str__
%ignore NumPyArrayBase::getDataArray() const;
%ignore NumPyArrayBase::getArrayObject() const;
%ignore Epetra_LinearProblem::SetOperator()(int) const;
%ignore Epetra_Time::operator=(const Epetra_Time &);

// Rename directives
%rename(Version             ) Epetra_Version;
%rename(Object              ) Epetra_Object;
%rename(BlockMap            ) Epetra_BlockMap;
%rename(Map                 ) Epetra_Map;
%rename(LocalMap            ) Epetra_LocalMap;
%rename(SrcDistObject       ) Epetra_SrcDistObject;
%rename(DistObject          ) Epetra_DistObject;
%rename(CompObject          ) Epetra_CompObject;
%rename(BLAS                ) Epetra_BLAS;
%rename(LAPACK              ) Epetra_LAPACK;
%rename(IntVector           ) Epetra_IntVector;
%rename(CrsGraph            ) Epetra_CrsGraph;
%rename(MapColoring         ) Epetra_MapColoring;
%rename(Operator            ) Epetra_Operator;
%rename(RowMatrix           ) Epetra_RowMatrix;
%rename(CrsMatrix           ) Epetra_CrsMatrix;
%rename(FECrsMatrix         ) Epetra_FECrsMatrix;
%rename(VbrMatrix           ) Epetra_VbrMatrix;
%rename(SerialDenseSolver   ) Epetra_SerialDenseSolver;
%rename(SerialDenseOperator ) Epetra_SerialDenseOperator;
%rename(SerialDenseMatrix   ) Epetra_SerialDenseMatrix;
%rename(SerialDenseVector   ) Epetra_SerialDenseVector;
%rename(IntSerialDenseMatrix) Epetra_IntSerialDenseMatrix;
%rename(IntSerialDenseVector) Epetra_IntSerialDenseVector;
%rename(SerialSymDenseMatrix) Epetra_SerialSymDenseMatrix;
%rename(LinearProblem       ) Epetra_LinearProblem;
%rename(Import              ) Epetra_Import;
%rename(Export              ) Epetra_Export;
%rename(Time                ) Epetra_Time;

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

%include "Epetra_Comm.i"

%include "Epetra_BlockMap.h"
%include "Epetra_Map.h"
%include "Epetra_LocalMap.h"
%include "Epetra_SrcDistObject.h"
%include "Epetra_DistObject.h"
%include "Epetra_CompObject.h"
%include "Epetra_BLAS.h"
%include "Epetra_LAPACK.h"

%include "Epetra_Vectors.i"

%include "Epetra_IntVector.h"
%include "Epetra_CrsGraph.h"
%include "Epetra_MapColoring.h"
%include "Epetra_Operator.h"
%include "Epetra_RowMatrix.h"
%include "Epetra_CrsMatrix.h"
%include "Epetra_FECrsMatrix.h"
%include "Epetra_VbrMatrix.h"
%include "Epetra_DataAccess.h"
%include "Epetra_SerialDenseSolver.h"
%include "Epetra_SerialDenseOperator.h"
%include "Epetra_SerialDenseMatrix.h"
%include "Epetra_SerialDenseVector.h"
%include "Epetra_SerialSymDenseMatrix.h"
%include "Epetra_IntSerialDenseMatrix.h"
%include "Epetra_IntSerialDenseVector.h"
%include "Epetra_LinearProblem.h"
%include "Epetra_CombineMode.h"
%include "Epetra_Import.h"
%include "Epetra_Export.h"
%include "Epetra_Time.h"

// Local interface includes
%include "NumPyArray.h"

// Extensions
%extend Epetra_Object {

  // Define the __str__() method, used by the python str() operator on any
  // object given to the python print command.
  string __str__() {
    stringstream os;
    self->Print(os);             // Put the output in os
    string s = os.str();         // Extract the string from os
    int last = s.length();       // Get the last index
    if (s.substr(last) == "\n")
      last-=1;                   // Ignore any trailing newline
    return s.substr(0,last);     // Return the string
  }

  // The Epetra_Object::Print(ostream) method is ignored and replaced by a
  // Print() method here that takes a python file as its argument.  If no
  // argument is given, then output is to standard out.
  PyObject * Print(PyObject*pf=NULL) const {
    if (pf == NULL) {
      self->Print(std::cout);
    } else {
      if (!PyFile_Check(pf)) {
	PyErr_SetString(PyExc_IOError, "Print() method expects file object");
	return NULL;
      } else {
	std::FILE*   f = PyFile_AsFile(pf);
	FileStream   buffer(f);
	std::ostream os(&buffer);
	self->Print(os);
      }
    }
    Py_INCREF(Py_None);
    return Py_None;
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
  void __setitem__(PyObject* args, double val) 
  {
    int Row, Col;
    if (!PyArg_ParseTuple(args, "ii", &Row, &Col)) {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return;
    }

    if (self->ReplaceGlobalValues(Row, 1, &val, &Col))
      self->InsertGlobalValues(Row, 1, &val, &Col);
  }

  PyObject* __getitem__(PyObject* args) 
  {
    int Row, Col;
    if (PyInt_Check(args))
    {
      return(Epetra_RowMatrix_GetEntries(*self, PyLong_AsLong(args)));
    }
    else if (PyArg_ParseTuple(args, "ii", &Row, &Col))
    {
      return(Epetra_RowMatrix_GetEntry(*self, Row, Col));
    }
    else
    {
      PyErr_SetString(PyExc_IndexError, "Input argument not supported");
      Py_INCREF(Py_None);
      return Py_None;
    }
  }

  int InsertGlobalValues(const int Row, const int Size, 
                         const Epetra_SerialDenseVector& Values,
                         const Epetra_IntSerialDenseVector& Entries)
  {
    return self->InsertGlobalValues(Row, Size, Values.Values(), (int*)Entries.Values());
  }

  int InsertGlobalValue(int i, int j, double val) {
    double val2 = val;
    int j2 = j;
    return self->InsertGlobalValues(i, 1, &val2, &j2);
  }


  int InsertGlobalValues(const int row, PyObject* Values, PyObject* Indices)
  {
    if (row < 0)
      return(-1);

    if (PyList_Check(Values) == 0 || PyList_Check(Indices) == 0) 
    {
      cerr << "Input object is not a list" << endl;
      return(-1);
    }

    int len = PyList_Size(Values);
    if (len != PyList_Size(Indices))
    {
      cerr << "Length of input lists differ" << endl;
      return(-1);
    }

    for (int i = 0 ; i < len ; ++i)
    {
      PyObject* Value,* Index;
      Value = PyList_GetItem(Values, i);
      Index = PyList_GetItem(Indices, i);

      if (PyInt_Check(Index) == 0)
      {
        cerr << "Indices must be integers" << endl;
        return(-1);
      }

      if (PyFloat_Check(Value) == 0)
      {
        cerr << "Values must be doubles" << endl;
        return(-1);
      }

      int cIndex = PyLong_AsLong(Index);
      double cValue = PyFloat_AsDouble(Value);
      if (self->InsertGlobalValues(row, 1, &cValue, &cIndex) < 0)
        return(-1);
    }
    return(0);
  }
}

%extend Epetra_FECrsMatrix {
  void __setitem__(PyObject* args, double val) 
  {
    int Row, Col;
    if (!PyArg_ParseTuple(args, "ii", &Row, &Col)) {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return;
    }

    if (self->ReplaceGlobalValues(1, &Row, 1, &Col, &val))
      self->InsertGlobalValues(1, &Row, 1, &Col, &val);
  }

  PyObject* __getitem__(PyObject* args) 
  {
    int Row, Col;
    if (PyInt_Check(args))
    {
      return(Epetra_RowMatrix_GetEntries(*self, PyLong_AsLong(args)));
    }
    else if (PyArg_ParseTuple(args, "ii", &Row, &Col))
    {
      return(Epetra_RowMatrix_GetEntry(*self, Row, Col));
    }
    else
    {
      PyErr_SetString(PyExc_IndexError, "Input argument not supported");
      Py_INCREF(Py_None);
      return Py_None;
    }
  }

  int InsertGlobalValues(const int Row, const int Size, 
                         const Epetra_SerialDenseVector& Values,
                         const Epetra_IntSerialDenseVector& Entries)
  {
    return self->InsertGlobalValues(1, &Row,
                                    Size, (int*)Entries.Values(),
                                    Values.Values());
  }

  int InsertGlobalValue(int i, int j, double val) {
    double val2 = val;
    int j2 = j;
    return self->InsertGlobalValues(1, &i, 1, &j2, &val2);
  }

  int InsertGlobalValues(const int row, PyObject* Values, PyObject* Indices)
  {
    if (row < 0)
      return(-1);

    if (PyList_Check(Values) == 0 || PyList_Check(Indices) == 0) 
    {
      cerr << "Input object is not a list" << endl;
      return(-1);
    }

    int len = PyList_Size(Values);
    if (len != PyList_Size(Indices))
    {
      cerr << "Length of input lists differ" << endl;
      return(-1);
    }

    for (int i = 0 ; i < len ; ++i)
    {
      PyObject* Value,* Index;
      Value = PyList_GetItem(Values, i);
      Index = PyList_GetItem(Indices, i);

      if (PyInt_Check(Index) == 0)
      {
        cerr << "Indices must be integers" << endl;
        return(-1);
      }

      if (PyFloat_Check(Value) == 0)
      {
        cerr << "Values must be doubles" << endl;
        return(-1);
      }

      int cIndex = PyLong_AsLong(Index);
      double cValue = PyFloat_AsDouble(Value);
      if (self->InsertGlobalValues(1, &row, 1, &cIndex, &cValue) < 0)
        return(-1);
    }
    return(0);
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

%extend Epetra_IntSerialDenseVector {

  int __call__(int i) {
    return self->operator()(i);
  }

  int __getitem__(int i) {
    return self->operator[](i);
  }

  void __setitem__(int i, const int val) {
    int * column = self->Values();
    column[i] = val;
  }
}

%extend Epetra_Map 
{
  Epetra_Map(const int NumGlobalElements,
             const Epetra_IntSerialDenseVector& MyGlobalElements,
             const int IndexBase, const Epetra_Comm& Comm)
  {
    return(new Epetra_Map(NumGlobalElements, MyGlobalElements.Length(),
                         (int*)MyGlobalElements.Values(), IndexBase, Comm));
  }

  Epetra_Map(const int NumGlobalElements,
             PyObject* MyGlobalElements, const int IndexBase,
             const Epetra_Comm& Comm)
  {
    if (PyList_Check(MyGlobalElements) == 0)
    {
      cerr << "Input object is not a list" << endl;
      return NULL;
    }

    int len = PyList_Size(MyGlobalElements);

    vector<int> list(len);

    for (int i = 0 ; i < len ; ++i)
    {
      PyObject* Index;
      Index = PyList_GetItem(MyGlobalElements, i);

      if (PyInt_Check(Index) == 0)
      {
        cerr << "Indices must be integers" << endl;
        return NULL;
      }

      list[i] = PyLong_AsLong(Index);
    }
    return(new Epetra_Map(NumGlobalElements, len, &list[0], IndexBase, Comm));
  }

  PyObject*  MyGlobalElements()
  {
    int* MyGlobalElements_Epetra = self->MyGlobalElements();
    PyObject* MyGlobalElements_Python,* item;
    int size = self->NumMyElements();
    if (size <= 0)
      goto fail;

    MyGlobalElements_Python = PyList_New(size);

    for (int i = 0 ; i < size ; ++i)
    {
      item = PyInt_FromLong(MyGlobalElements_Epetra[i]);
      PyList_SetItem(MyGlobalElements_Python, i, item);
    }

    return(MyGlobalElements_Python);
fail:
    Py_INCREF(Py_None);
    return Py_None;
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

// Python code.  Here we set the __version__ string, call MPI_Init()
// and arrange for MPI_Finalize to be called (if appropriate), and
// declare classes that inherit both from Epetra objects and
// UserArrays, to give users additional functionality.
%pythoncode %{

__version__ = Version().split()[2]

# Call MPI_Init if appropriate
import sys
Init_Argv(sys.argv)
del sys

# Arrange for MPI_Finalize to be called at exit, if appropriate
import atexit
atexit.register(Finalize)

from UserArray import *

class MultiVector(UserArray,NumPyMultiVector):
    def __init__(self, *args):
        """
        __init__(self, BlockMap map, int numVectors, bool zeroOut=True) -> MultiVector
        __init__(self, MultiVector source) -> MultiVector
        __init__(self, BlockMap map, PyObject array) -> MultiVector
        __init__(self, DataAccess CV, MultiVector source, PyObject range) -> MultiVector
        __init__(self, PyObject array) -> MultiVector
        """
        NumPyMultiVector.__init__(self, *args)
        UserArray.__init__(self,self.getArray(),'d',copy=False,savespace=False)
    def __str__(self):
        return str(self.array)
    def __setattr__(self, key, value):
        "Protect the 'array' and 'shape' attributes"
        if key == "array":
            if key in self.__dict__:
                raise AttributeError, "Cannot change Epetra.MultiVector array attribute"
        elif key == "shape":
            value = tuple(value)
            if len(value) < 2:
                raise ValueError, "Epetra.MultiVector shape is " + str(value) + \
		  " but must have minimum of 2 elements"
        UserArray.__setattr__(self, key, value)

class Vector(UserArray,NumPyVector):
    def __init__(self, *args):
        NumPyVector.__init__(self, *args)
        UserArray.__init__(self,self.getArray(),'d',copy=False,savespace=False)
    def __str__(self):
        return str(self.array)
    def __setattr__(self, key, value):
        "Protect the 'array' attribute"
        if key == "array":
            if key in self.__dict__:
                raise AttributeError, "Cannot change Epetra.Vector array attribute"
        UserArray.__setattr__(self, key, value)

%}

#ifndef HAVE_MPI
%pythoncode %{
def PyComm():
  return SerialComm();
%}
#else
%pythoncode %{
def PyComm():
  return MpiComm(CommWorld());
%}
#endif

%include "Epetra_PyOperator.h"
%include "Epetra_PyRowMatrix.h"
