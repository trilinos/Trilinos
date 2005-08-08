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

%define ML_DOCSTRING
"The ML module allows access to The Trilinos package ML.  Note
that the 'ML_' prefix has been stripped from all ML objects,
but that if imported with 'from PyTrilinos import ML', these
objects exist in the 'ML' python namespace.  Use the python help()
facility for local documentation on classes and methods, or see the
on-line documentation for more in-depth information.

The most important class of the ML module is:
- MultiLevelPreconditioner

Example of usage:
1) Creates an ML preconditioner for the already created Epetra.CrsMatrix 
   Matrix. The parameters are set using a dictionary. 

   Prec = ML.MultiLevelPreconditioner(Matrix, False)
   MLList = {
     \"max levels\"        : 3,
     \"smoother: type\"    : \"symmetric Gauss-Seidel\",
     \"aggregation: type\" : \"Uncoupled\"
   }
   Prec.SetParameterList(MLList)
   Prec.ComputePreconditioner()
"
%enddef

%module(package="PyTrilinos", directors="1", docstring=ML_DOCSTRING) ML

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

#include "ml_config.h"

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_PyOperator.h"
#include "Epetra_PyRowMatrix.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_ParameterList.hpp"
#include "MLAPI.h"
#include "MLAPI_PyMatrix.h"

Teuchos::ParameterList* CreateList(PyObject* obj)
{
  int i;
  Teuchos::ParameterList* List;
  if (!PyDict_Check(obj)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a dictionary");
    return NULL;
  }
  List = new Teuchos::ParameterList;

  int size = PyDict_Size(obj);
  PyObject* Keys = PyDict_Keys(obj);
  PyObject* Values = PyDict_Values(obj);

  for (i = 0; i < size ; i++) 
  {
    PyObject *s = PyList_GetItem(Keys,i);
    PyObject *t = PyList_GetItem(Values,i);

    // Get the parameter name
    if (!PyString_Check(s)) {
        PyErr_SetString(PyExc_ValueError, "Dictionary keys must be strings");
        return NULL;
    }
    string ParameterName = PyString_AsString(s);

    // now parse for the parameter value and type
    // This can be a "int", "double", "string", or a tuple
    // for more general types

    if (PyBool_Check(t)) 
    {
      if (t == Py_True)
        List->set(ParameterName, true);
      else
        List->set(ParameterName, false);
    }
    else if (PyInt_Check(t)) 
    {
      int ParameterValue = PyInt_AsLong(t);
      List->set(ParameterName, ParameterValue);
    }
    else if (PyFloat_Check(t)) 
    {
      double ParameterValue = PyFloat_AsDouble(t);
      List->set(ParameterName, ParameterValue);
    }
    else if (PyString_Check(t)) 
    {
      string ParameterValue = PyString_AsString(t);
      List->set(ParameterName, ParameterValue);
    }
    else if (PyTuple_Check(t)) 
    {
      if (!PyString_Check(PyTuple_GetItem(t, 0)) ||
          !PyString_Check(PyTuple_GetItem(t, 1))) {
        PyErr_SetString(PyExc_ValueError, "tuples must contain strings");
        return NULL;
      }
      string ParameterType = PyString_AsString(PyTuple_GetItem(t, 0));
      string ParameterValue = PyString_AsString(PyTuple_GetItem(t, 1));
      if (ParameterType == "bool") 
      {
        if (ParameterValue == "true")
          List->set(ParameterName, true);
        else
          List->set(ParameterName, false);
      }
      else if (ParameterType == "int") 
      {
        List->set(ParameterName, (int)atoi(ParameterValue.c_str()));
      }
      else if (ParameterType == "double") 
      {
        List->set(ParameterName, (double)atof(ParameterValue.c_str()));
      }
      else if (ParameterType == "string") 
      {
        List->set(ParameterName, string(ParameterValue));
      }
      else 
      {
        PyErr_SetString(PyExc_ValueError, "type in tuple not recognized");
        return NULL;
      }
    }
    else
    {
      PyErr_SetString(PyExc_ValueError, "Type in list not recognized");
      return NULL;
    }
  }

  return(List);
}

PyObject * Epetra_RowMatrix_GetEntries(const Epetra_RowMatrix& Matrix,
                                       int GlobalRow) 
{
  if (!Matrix.Filled())
  {
    cout << "Matrix not FillComplete()'d" << endl;
    Py_INCREF(Py_None);
    return Py_None;
  }
  int ierr;
  PyObject* res;
  PyObject* PyIndices;
  PyObject* PyValues;
  int NumEntries, Length;
  Length = Matrix.MaxNumEntries();
  vector<int>    Indices(Length);
  vector<double> Values(Length);
  int MyRow = Matrix.RowMatrixRowMap().LID(GlobalRow);
  ierr = Matrix.ExtractMyRowCopy(MyRow, Length, NumEntries, &Values[0],
                                 &Indices[0]);
  if (ierr < 0)
  {
    Py_INCREF(Py_None);
    return Py_None;
  }

  PyIndices = PyList_New(NumEntries);
  PyValues  = PyList_New(NumEntries);

  // return global indices
  for (int i = 0 ; i < NumEntries ; ++i)
  {
    int GID = Matrix.RowMatrixColMap().GID(Indices[i]);
    PyList_SetItem(PyIndices, i, PyInt_FromLong(GID));
    PyList_SetItem(PyValues,  i, PyFloat_FromDouble(Values[i]));
  }
  res = PyTuple_New(2);
  PyTuple_SetItem(res, 0, PyIndices);
  PyTuple_SetItem(res, 1, PyValues);
  return(res);
}

PyObject* Epetra_RowMatrix_GetEntry(Epetra_RowMatrix& Matrix, 
                                    int GlobalRow, int GlobalCol)
{
  int NumEntries, Length, ierr;
  double val = 0.0;
  Length = Matrix.MaxNumEntries();
  vector<int>    Indices(Length);
  vector<double> Values(Length);
  int MyRow = Matrix.RowMatrixRowMap().LID(GlobalRow);
  int MyCol = Matrix.RowMatrixColMap().LID(GlobalCol);

  ierr = Matrix.ExtractMyRowCopy(MyRow, Length, NumEntries, &Values[0], 
                                 &Indices[0]);
  if (ierr) NumEntries = 0;
  for (int i = 0 ; i < Length ; ++i)
  {
    if (Indices[i] == MyCol)
    {
      val = Values[i];
      break;
    }
  }
  return(PyFloat_FromDouble(val));
}

MLAPI::Operator GetPNonSmoothed(const MLAPI::Operator& A,
                                const MLAPI::MultiVector& ThisNS,
                                MLAPI::MultiVector& NextNS,
                                PyObject* obj)
{   
  Teuchos::ParameterList* List = CreateList(obj);
  MLAPI::Operator Ptent;
  MLAPI::GetPtent(A, *List, ThisNS, Ptent, NextNS);
  delete List;
  return(Ptent);
}

bool Iterate(const MLAPI::Operator& A, const MLAPI::MultiVector& LHS,
             const MLAPI::MultiVector& RHS, const MLAPI::BaseOperator& Prec, 
             PyObject* obj)
{
  Teuchos::ParameterList* List;
  List = CreateList(obj);
  if (List == 0)
    return(false);
  Krylov(A, LHS, RHS, Prec, *List);
  delete List;
  return(true);
}

extern "C" {
  // on some MAC OS X with LAM/MPI _environ() is not found,
  // need to specify -Wl,-i_environ:_fake_environ as LDFLAGS
  void fake_environ()
  {
    exit(EXIT_FAILURE);
  }
}

%}

using namespace std;
%include "ml_config.h"

// Auto-documentation feature
%feature("autodoc", "1");

%typemap(in) (const ParameterList & List)
{
  $1 = CreateList($input);
  cout << *($1);
}

%typemap(in) (Teuchos::ParameterList& List)
{
  $1 = CreateList($input);
  cout << *($1);
}

%typemap(freearg) (const MLAPI::Operator& A, Teuchos::ParameterList& List,
                   const MLAPI::MultiVector& ThisNS, 
                   MLAPI::Operator& Ptent, MLAPI::MultiVector& NextNS)
{   
  if ($2) delete($2);
}   

%feature("director") MLAPI::BaseOperator;
%import "Epetra.i"

// Amesos interface includes
%include "ml_MultiLevelPreconditioner.h"

%include "MLAPI_Workspace.h"
%include "MLAPI_BaseObject.h"
%include "MLAPI_CompObject.h"
%include "MLAPI_TimeObject.h"
%include "MLAPI_BaseOperator.h"
%include "MLAPI_Space.h"
%include "MLAPI_MultiVector.h"
%include "MLAPI_MultiVector_Utils.h"
%include "MLAPI_Operator.h"
%include "MLAPI_InverseOperator.h"
%include "MLAPI_Operator_Utils.h"
%include "MLAPI_EpetraBaseOperator.h"
%include "MLAPI_Krylov.h"
%include "MLAPI_Expressions.h"
%include "MLAPI_Aggregation.h"
%include "MLAPI_InverseOperator.h"
%include "MLAPI_Eig.h"
%include "MLAPI_Gallery.h"
%include "MLAPI_PyMatrix.h"

%extend MLAPI::Space {

  PyObject*  GetMyGlobalElements()
  {
    int n = self->GetNumMyElements();
    if (n <= 0)
      goto fail;

    PyObject* MyGlobalElements_Python,* item;

    MyGlobalElements_Python = PyList_New(n);

    for (int i = 0 ; i < n ; ++i)
    {
      item = PyInt_FromLong((*self)(i));
      PyList_SetItem(MyGlobalElements_Python, i, item);
    }

    return(MyGlobalElements_Python);
fail:
    Py_INCREF(Py_None);
    return Py_None;
  }
}

%extend MLAPI::BaseObject {

  // Define the __str__() method, used by the python str() operator on any
  // object given to the python print command.
  string __str__() {
    stringstream os;
    self->Print(os);                  // Put the output in os
    string s = os.str();              // Extract the string from os
    return s.substr(0,s.length()-1);  // Return the string minus trailing \n
  }
}

%extend MLAPI::Operator {

  PyObject* __getitem__(PyObject* args) 
  {
    Epetra_RowMatrix* Matrix = self->GetRCPRowMatrix().get();

    int Row, Col;
    if (PyInt_Check(args))
    {
      return(Epetra_RowMatrix_GetEntries(*Matrix, PyLong_AsLong(args)));
    }
    else if (PyArg_ParseTuple(args, "ii", &Row, &Col))
    {
      return(Epetra_RowMatrix_GetEntry(*Matrix, Row, Col));
    }
    else
    {
      PyErr_SetString(PyExc_IndexError, "Input argument not supported");
      Py_INCREF(Py_None);
      return Py_None;
    }
  }

  MLAPI::MultiVector __mul__(MLAPI::MultiVector& rhs)
  {
    return(*self * rhs);
  }

  MLAPI::Operator __add__(MLAPI::Operator& rhs)
  {
    return(*self + rhs);
  }

  MLAPI::Operator __sub__(MLAPI::Operator& rhs)
  {
    return(*self - rhs);
  }

  MLAPI::Operator __mul__(MLAPI::Operator& rhs)
  {
    // questo non lo capisco ma senza non funziona...
    if (self == &rhs)
      return(*self * MLAPI::Duplicate(rhs));
    else
      return(*self * rhs);
  }

  MLAPI::Operator __mul__(double rhs)
  {
    return(*self * rhs);
  }

  MLAPI::Operator __div__(double rhs)
  {
    return(*self * (1.0 / rhs));
  }
}

%extend MLAPI::InverseOperator {
  MLAPI::MultiVector __mul__(MLAPI::MultiVector& rhs)
  {
    return(*self * rhs);
  }

  bool Reshape(const Operator& Op, const string Type, PyObject* obj)
  {
    Teuchos::ParameterList* List = CreateList(obj);
    if (List == 0)
      return(false);
    else
      self->Reshape(Op, Type, *List);
    delete List;
    return(true);
  }
}

%extend MLAPI::MultiVector {
  MLAPI::MultiVector __add__(MLAPI::MultiVector& rhs)
  {
    return(*self + rhs);
  }

  MLAPI::MultiVector __sub__(MLAPI::MultiVector& rhs)
  {
    return(*self - rhs);
  }

  double __mul__(MLAPI::MultiVector& rhs)
  {
    return(*self * rhs);
  }

  MLAPI::MultiVector __mul__(double rhs)
  {
    return(*self * rhs);
  }

  void __setitem__(PyObject* args, double value)
  {
    int Row, Col;
    if (!PyArg_ParseTuple(args, "ii", &Row, &Col)) {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return;
    }

    (*self)(Row, Col) = value;
  }

  double __getitem__(PyObject* args)
  {
    int Row, Col;
    if (!PyArg_ParseTuple(args, "ii", &Row, &Col)) {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return(0.0);
    }

    return((*self)(Row, Col));
  }
}

%extend PyMatrix {
  void __setitem__(PyObject* args, double val) 
  {
    int Row, Col;
    if (!PyArg_ParseTuple(args, "ii", &Row, &Col)) {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return;
    }

    self->SetElement(Row, Col, val);
  }
  
  PyObject* __getitem__(PyObject* args) 
  {
    int Row, Col;
    if (PyInt_Check(args))
    {
      return(Epetra_RowMatrix_GetEntries(*(self->GetMatrix()), PyLong_AsLong(args)));
    }
    else if (PyArg_ParseTuple(args, "ii", &Row, &Col))
    {
      return(Epetra_RowMatrix_GetEntry(*(self->GetMatrix()), Row, Col));
    }
    else
    {
      PyErr_SetString(PyExc_IndexError, "Input argument not supported");
      Py_INCREF(Py_None);
      return Py_None;
    }
  }
}

MLAPI::Operator GetPNonSmoothed(const MLAPI::Operator& A,
                                const MLAPI::MultiVector& ThisNS,
                                MLAPI::MultiVector& NextNS,
                                PyObject* obj);
bool Iterate(const MLAPI::Operator& A, const MLAPI::MultiVector& LHS,
             const MLAPI::MultiVector& RHS, const MLAPI::BaseOperator& Prec, 
             PyObject* obj);

%pythoncode %{

# Call MPI_Init if appropriate
Init()

# Arrange for MPI_Finalize to be called at exit, if appropriate
import atexit
atexit.register(Finalize)

%}
