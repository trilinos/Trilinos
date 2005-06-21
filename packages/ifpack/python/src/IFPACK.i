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

%module(package="PyTrilinos") IFPACK

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

#include "Ifpack_ConfigDefs.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_VbrMatrix.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Operator.h"
#include "Teuchos_ParameterList.hpp"

// Amesos includes
#include "Ifpack.h"
#include "Ifpack_Utils.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_IC.h"
#include "Ifpack_ICT.h"
#include "Ifpack_ILU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_PointRelaxation.h"

extern "C" {
  // on some MAC OS X with LAM/MPI _environ() is not found,
  // need to specify -Wl,-i_environ:_fake_environ as LDFLAGS
  void fake_environ()
  {
    exit(EXIT_FAILURE);
  }
}

%}

// Auto-documentation feature
%feature("autodoc", "1");

%rename(Factory               ) Ifpack;
%rename(Preconditioner        ) Ifpack_Preconditioner;
%rename(IC                    ) Ifpack_IC;
%rename(ICT                   ) Ifpack_ICT;
%rename(ILU                   ) Ifpack_ILU;
%rename(ILUT                  ) Ifpack_ILUT;
%rename(Amesos                ) Ifpack_Amesos;
%rename(PointRelaxation       ) Ifpack_PointRelaxation;
%rename(AnalyzeMatrix         ) Ifpack_Analyze;
%rename(AnalyzeMatrixElements ) Ifpack_AnalyzeMatrixElements;
%rename(AnalyzeVectorElements ) Ifpack_AnalyzeVectorElements;
%rename(PrintSparsity         ) Ifpack_PrintSparsity;

%include "Ifpack_config.h"
%include "Ifpack_ConfigDefs.h"

// typemaps
%typemap(in) (Teuchos::ParameterList& List)
{
  int i;
  if (!PyDict_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a dictionary");
    return NULL;
  }
  $1 = new Teuchos::ParameterList;

  int size = PyDict_Size($input);
  PyObject* Keys = PyDict_Keys($input);
  PyObject* Values = PyDict_Values($input);

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
        $1->set(ParameterName, true);
      else
        $1->set(ParameterName, false);
    }
    else if (PyInt_Check(t)) 
    {
      int ParameterValue = PyInt_AsLong(t);
      $1->set(ParameterName, ParameterValue);
    }
    else if (PyFloat_Check(t)) 
    {
      double ParameterValue = PyFloat_AsDouble(t);
      $1->set(ParameterName, ParameterValue);
    }
    else if (PyString_Check(t)) 
    {
      string ParameterValue = PyString_AsString(t);
      $1->set(ParameterName, ParameterValue);
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
          $1->set(ParameterName, true);
        else
          $1->set(ParameterName, false);
      }
      else if (ParameterType == "int") 
      {
        $1->set(ParameterName, (int)atoi(ParameterValue.c_str()));
      }
      else if (ParameterType == "double") 
      {
        $1->set(ParameterName, (double)atof(ParameterValue.c_str()));
      }
      else if (ParameterType == "string") 
      {
        $1->set(ParameterName, string(ParameterValue));
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
}

%typemap(freearg) Teuchos::ParameterList& List
{
  delete($1);
}

using namespace std;

// Epetra interface includes
%import "Epetra.i"

// IFPACK interface includes
%include "Ifpack.h"
%include "Ifpack_Utils.h"
%include "Ifpack_Preconditioner.h"
%include "Ifpack_IC.h"
%include "Ifpack_ICT.h"
%include "Ifpack_ILU.h"
%include "Ifpack_ILUT.h"
%include "Ifpack_PointRelaxation.h"

%extend Ifpack_Preconditioner
{
  string __str__() {
    stringstream os;
    os << *self;
    return os.str();
  }

  void __del__()
  {
    delete self;
  }
}

