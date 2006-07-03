// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//                PyTrilinos.ML: Python Interface to ML
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
"The ML module allows access to The Trilinos package ML/MLAPI. Note
that the 'ML_' prefix has been stripped from all ML objects,
but that if imported with 'from PyTrilinos import ML', these
objects exist in the 'ML' python namespace.  Use the python help()
facility for local documentation on classes and methods, or see the
on-line documentation for more in-depth information.


Brief Description
=================

The ML module offers most of ML and MLAPI through Python.

The most important classes of the ML/MLAPI module are:
- MultiLevelPreconditioner
- PyMatrix
- MultiVector
- BaseOperator
- Operator
- InverseOperator
- EpetraBaseOperator
- <Operator> + <Operator>
- <Operator> * <Operator>
- <Operator> * <MultiVector>
- GetIdentity, GetTranspose, Iterate, ...


Example of usage: The MultiLevelPreconditioner class
====================================================

Creates an ML preconditioner for the already created Epetra.CrsMatrix 
Matrix. The parameters are set using a dictionary. 

from PyTrilinos Epetra, ML
<create here Matrix>
Prec = ML.MultiLevelPreconditioner(Matrix, False)
MLList = {
  \"max levels\"        : 3,
  \"smoother: type\"    : \"symmetric Gauss-Seidel\",
  \"aggregation: type\" : \"Uncoupled\"
}
Prec.SetParameterList(MLList)
Prec.ComputePreconditioner()


Example of usage: The MLAPI Interface
=====================================

from PyTrilinos import Epetra, ML
n = 1000
Space = ML.Space(n)
        
A = ML.PyMatrix(Space, Space)
      
for i in Space.GetMyGlobalElements():
  if i > 0:
    A[i, i - 1] = -1.
  if i < n - 1:
    A[i, i + 1] = -1.
  A[i, i] = 2.0
                               
A.FillComplete()

x = ML.MultiVector(Space); x.Random()
y = ML.MultiVector(Space)
r = y - A * x
print r.Norm2()
"
%enddef

%module(package="PyTrilinos", directors="1", docstring=ML_DOCSTRING) ML

// This is to avoid BaseLinearCombination and derived classes.
// MLAPI_LC is defined in setup.py; the MLAPI code contains some 
// `#ifndef MLAPI_LC' that must be kept in place.
#define MLAPI_LC

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

#include "ml_config.h"

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_FEVector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_JadOperator.h"

// Epetra python includes
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_PyOperator.h"
#include "Epetra_PyRowMatrix.h"

#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_ParameterList.hpp"
#include "MLAPI.h"
#include "MLAPI_PyMatrix.h"
#include "PyEpetra_Utils.h"
#include "PyTeuchos_Utils.h"

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

%}

// Ignore directives
%ignore *::operator=;
%ignore *::operator[];

using namespace std;
%include "ml_config.h"

// Auto-documentation feature
%feature("autodoc", "1");

%typemap(in) (const ParameterList & List)
{
  $1 = CreateList($input);
}

%typemap(in) (Teuchos::ParameterList& List)
{
  $1 = CreateList($input);
}

%typemap(freearg) (const MLAPI::Operator& A, Teuchos::ParameterList& List,
                   const MLAPI::MultiVector& ThisNS, 
                   MLAPI::Operator& Ptent, MLAPI::MultiVector& NextNS)
{   
  if ($2) delete($2);
}   

%feature("director") MLAPI::BaseOperator;
%import "Epetra.i"

// ML interface includes
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

%extend ML_Epetra::MultiLevelPreconditioner
{
  // This is for MatrixPortal
  int SetParameterListAndNullSpace(PyObject* obj,
                                   Epetra_MultiVector& NullSpace)
  {
    Teuchos::ParameterList* List = CreateList(obj);

    // WARNING: THIS IS DELICATE, NULLSPACE SHOULD NOT DISAPPEAR
    // otherwise the pointer here stored will vanish. This function should
    // be used only through the MatrixPortal
    double* NullSpacePtr = (double*)NullSpace.Values();
    int NullSpaceDim = NullSpace.NumVectors();

    List->set("null space: type", "pre-computed");
    List->set("null space: vectors", NullSpacePtr);
    List->set("null space: dimension", NullSpaceDim);

    self->SetParameterList(*List);

    delete List;

    return(0);
  }
}

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
    MultiVector res = *self + rhs;

    return(res);
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
