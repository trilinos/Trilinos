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

%{
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_CrsSingletonFilter.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_JadOperator.h"
#include "Epetra_LinearProblem.h"

// Local includes
#include "Epetra_PyOperator.h"
#include "Epetra_PyRowMatrix.h"
%}

// Feature directives
%feature("director") Epetra_Operator;
%feature("director") PyOperator;
%feature("director") Epetra_RowMatrix;
%feature("director") PyRowMatrix;

// Ignore directives
%ignore Epetra_LinearProblem::SetOperator()(int) const;
%ignore Epetra_VbrMatrix::Solve(bool, bool, bool,
				Epetra_Vector const&, Epetra_Vector&) const;

// Rename directives
%rename(Operator          ) Epetra_Operator;
%rename(RowMatrix         ) Epetra_RowMatrix;
%rename(CrsMatrix         ) Epetra_CrsMatrix;
%rename(FECrsMatrix       ) Epetra_FECrsMatrix;
%rename(CrsSingletonFilter) Epetra_CrsSingletonFilter;
%rename(VbrMatrix         ) Epetra_VbrMatrix;
%rename(FEVbrMatrix       ) Epetra_FEVbrMatrix;
%rename(JadOperator       ) Epetra_JadOperator;
%rename(LinearProblem     ) Epetra_LinearProblem;

// Exceptions
EXCEPTION_HANDLER(Epetra_CrsMatrix    ,OptimizeStorage   )
EXCEPTION_HANDLER(Epetra_FastCrsMatrix,FastCrsMatrix     )
EXCEPTION_HANDLER(Epetra_JadOperator  ,Epetra_JadOperator)

// Typemap directives

// // Begin input typemap collection for (int NumEntries, double * Values)
// %typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (int NumEntries, double * Values) {
//   $1 = ($input != 0);
// }
// %typemap(in) (int NumEntries, double * Values) %{
//   int entryLengths = 0;
//   PyArrayObject * valArray = (PyArrayObject*) PyArray_ContiguousFromObject($input,
// 									   'd', 0, 0);
//   if (valArray == NULL) SWIG_exception(SWIG_ValueError,"Invalid sequence of values");
//   entryLengths = _PyArray_multiply_list(valArray->dimensions,valArray->nd);
//   $1 = entryLengths;
//   $2 = (double *) (valArray->data);
// %}
// %typemap(freearg) (int NumEntries, double * Values) %{
//   Py_XDECREF(valArray);
// %}
// // End input typemap collection for (int NumEntries, double * Values)

// // Begin input typemap collection for (int * Indices)
// %typecheck(SWIG_TYPECHECK_INT32_ARRAY) (int * Indices) {
//   $1 = ($input != 0);
// }
// %typemap(in) (int * Indices) %{
//   PyArrayObject * indArray = (PyArrayObject*) PyArray_ContiguousFromObject($input,
// 									   'i',0,0);
//   if (indArray == NULL) SWIG_exception(SWIG_ValueError,"Invalid sequence of indices");
//   if (_PyArray_multiply_list(indArray-dimensions,indArray->nd) != entryLengths)
//     SWIG_exception(SWIG_ValueError,"Values and Indices have different lengths");
//   $1 = (int*) (indArray->data);
// %}
// %typemap(freearg) (int * Indices) %{
//   Py_XDECREF(indArray);
// %}
// // End input typemap collection for (int * Indices)

// Include directives
%warnfilter(473) Epetra_Operator;
%warnfilter(473) Epetra_RowMatrix;
%include "Epetra_Operator.h"
%include "Epetra_RowMatrix.h"
%include "Epetra_CrsMatrix.h"
%include "Epetra_FECrsMatrix.h"
%include "Epetra_CrsSingletonFilter.h"
%include "Epetra_VbrMatrix.h"
%include "Epetra_FEVbrMatrix.h"
%include "Epetra_JadOperator.h"
%include "Epetra_LinearProblem.h"

// Clear the typemaps
%clear (int NumEntries, int * Values);
%clear (int * Indices);

// Extend directives
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

%warnfilter(473) PyOperator;
%warnfilter(473) PyRowMatrix;
%include "Epetra_PyOperator.h"
%include "Epetra_PyRowMatrix.h"

// Epetra_RowMatrixTransposer, Epetra_FastCrsMatrix and
// Epetra_LinearProblemRedistor are apparently not built
//#include "Epetra_RowMatrixTransposer.h"
//#include "Epetra_FastCrsMatrix.h"
//#include "Epetra_LinearProblemRedistor.h"
//%rename(RowMatrixTransposer  ) Epetra_RowMatrixTransposer;
//%rename(FastCrsMatrix        ) Epetra_FastCrsMatrix;
//%rename(LinearProblemRedistor) Epetra_LinearProblem;
//%include "Epetra_RowMatrixTransposer.h"
//%include "Epetra_FastCrsMatrix.h"
//%include "Epetra_LinearProblemRedistor.h"
