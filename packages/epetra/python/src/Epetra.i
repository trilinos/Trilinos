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

An Epetra.Vector object is at the same time a Numeric array, and it
can therefore be used everywhere Numeric arrays are accepted.
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
#include "Epetra_SrcDistObject.h"
#include "Epetra_DistObject.h"
#include "Epetra_CompObject.h"
#include "Epetra_BLAS.h"
#include "Epetra_LAPACK.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_MapColoring.h"
#include "Epetra_DataAccess.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
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

// tools for handling double* and int* pointers in Python
#include "Epetra_IArray.h"
#include "Epetra_DArray.h"
%}

%feature("director") PyOperator;
%feature("director") PyRowMatrix;

// Ignore directives
%ignore *::operator=;         // Not overrideable in python
%ignore *::operator[];        // Replaced with __setitem__ method
%ignore *::operator[] const;  // Replaced with __getitem__ method
%ignore *::print;             // Replaced with __str__ method
%ignore operator<<(ostream &, const Epetra_Object &);// From python, use __str__
%ignore NumPyArrayBase::getDataArray() const;
%ignore NumPyArrayBase::getArrayObject() const;
%ignore Epetra_Object::Print(ostream &) const;       // Replaced with __str__ method
%ignore Epetra_MapColoring::operator()(int) const;
%ignore Epetra_CompObject::UpdateFlops(int) const;   // Use long int version
%ignore Epetra_CompObject::UpdateFlops(float) const; // Use double version
%ignore Epetra_LinearProblem::SetOperator()(int) const;
%ignore Epetra_VbrMatrix::Solve(bool, bool, bool,
				Epetra_Vector const&, Epetra_Vector&) const;
// Epetra_VbrMatrix member function
// Solve(bool,bool,bool,Epetra_Vector const&,Epetra_Vector&) const does not
// appear to be implemented.  Apparently overridden by
// Solve(bool,bool,bool,Epetra_MultiVector const&,Epetra_MultiVector&) const

// Rename directives
%rename(Version      ) Epetra_Version;
%rename(Object       ) Epetra_Object;
%rename(SrcDistObject) Epetra_SrcDistObject;
%rename(DistObject   ) Epetra_DistObject;
%rename(CompObject   ) Epetra_CompObject;
%rename(BLAS         ) Epetra_BLAS;
%rename(LAPACK       ) Epetra_LAPACK;
%rename(CrsGraph     ) Epetra_CrsGraph;
%rename(MapColoring  ) Epetra_MapColoring;
%rename(Operator     ) Epetra_Operator;
%rename(RowMatrix    ) Epetra_RowMatrix;
%rename(CrsMatrix    ) Epetra_CrsMatrix;
%rename(FECrsMatrix  ) Epetra_FECrsMatrix;
%rename(VbrMatrix    ) Epetra_VbrMatrix;
%rename(LinearProblem) Epetra_LinearProblem;
%rename(Import       ) Epetra_Import;
%rename(Export       ) Epetra_Export;
%rename(Time         ) Epetra_Time;

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

// Epetra interface includes and imports
using namespace std;
%include "Epetra_Version.h"
%include "Epetra_CombineMode.h"
%include "Epetra_DataAccess.h"
%include "Epetra_Object.h"
%include "Epetra_SrcDistObject.h"
%include "Epetra_DistObject.h"
%include "Epetra_CompObject.h"
%import  "Epetra_BLAS.h"       // These two classes are not included because I do not
%import  "Epetra_LAPACK.h"     // want to expose their functionality to python
%include "Epetra_Time.h"
%include "Epetra_Import.h"
%include "Epetra_Export.h"
%include "Epetra_MapColoring.h"

%include "Epetra_SerialDense.i"
%include "Epetra_Comm.i"
%include "Epetra_Maps.i"
%include "Epetra_Vectors.i"
//%include Epetra_Graphs.i
//%include Epetra_Operators.i

%include "Epetra_CrsGraph.h"
%include "Epetra_Operator.h"
%include "Epetra_RowMatrix.h"
%include "Epetra_CrsMatrix.h"
%include "Epetra_FECrsMatrix.h"
%include "Epetra_VbrMatrix.h"
%include "Epetra_LinearProblem.h"

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

%extend NumPyArrayBase {
  using namespace std;
  string __str__() {
    stringstream os;
    self->print(os);                  // Put the output in os
    string s = os.str();              // Extract the string from os
    return s.substr(0,s.length()-1);  // Return the string minus trailing \n
  }
}

// Python code.  Here we set the __version__ string
%pythoncode %{

__version__ = Version().split()[2]

%}

%warnfilter(473) PyOperator;
%warnfilter(473) PyRowMatrix;
%include "Epetra_PyOperator.h"
%include "Epetra_PyRowMatrix.h"

%include "Epetra_IArray.h"
%include "Epetra_DArray.h"

%extend IArray {
  int __getitem__(int i)
  {
    return((*self)[i]);
  }

  void __setitem__(int i, int val)
  {
    (*self)[i] = val;
  }
}

%extend DArray {
  double __getitem__(int i)
  {
    return((*self)[i]);
  }

  void __setitem__(int i, double val)
  {
    (*self)[i] = val;
  }
}
