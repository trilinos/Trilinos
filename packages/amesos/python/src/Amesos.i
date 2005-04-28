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

%module(package="PyTrilinos") Amesos

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

#include "Amesos_ConfigDefs.h"

// Epetra includes
#include "Epetra_Object.h"
#include "Epetra_DistObject.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_DataAccess.h"

// Local includes
#include "Epetra_NumPyVector.h"
#include "NumPyArray.h"
#include "NumPyWrapper.h"

// Teuchos includes
#include "Teuchos_ParameterList.hpp"

// Amesos includes
#include "Amesos_ConfigDefs.h"
#include "Amesos.h"
#include "Amesos_BaseSolver.h"
#ifdef HAVE_AMESOS_LAPACK
#include "Amesos_Lapack.h"
#endif
#ifdef HAVE_AMESOS_KLU
#include "Amesos_Klu.h"
#endif
#ifdef HAVE_AMESOS_UMFPACK
#include "Amesos_Umfpack.h"
#endif
#ifdef HAVE_AMESOS_SUPERLU
#include "Amesos_Superlu.h"
#endif
#ifdef HAVE_AMESOS_MUMPS
#include "Amesos_Mumps.h"
#endif
%}

// Auto-documentation feature
%feature("autodoc", "1");

%include "Amesos_config.h"
%include "Amesos_ConfigDefs.h"

// Epetra
%import "Epetra.i"

// Rename directives for Amesos
%rename(BaseSolver          ) Amesos_BaseSolver;
%rename(Factory             ) Amesos;
%rename(Klu                 ) Amesos_Klu;
%rename(Lapack              ) Amesos_Lapack;
%rename(Umfpack             ) Amesos_Umfpack;
%rename(Superlu             ) Amesos_Superlu;
%rename(Superludist         ) Amesos_Superludist;
%rename(Mumps               ) Amesos_Mumps;
%rename(Dscpack             ) Amesos_Dscpack;

// SWIG library includes
%include "std_string.i"

// typemaps
%typemap(in) (Teuchos::ParameterList& ParameterList)
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

  for (i = 0; i < size ; i++) {
    PyObject *s = PyList_GetItem(Keys,i);
    PyObject *t = PyList_GetItem(Values,i);
    if (!PyString_Check(s)) {
        PyErr_SetString(PyExc_ValueError, "Dictionary keys must be strings");
        return NULL;
    }
    if (!PyTuple_Check(t)) {
        PyErr_SetString(PyExc_ValueError, "Dictionary values must be tuples");
        return NULL;
    }
    if (!PyString_Check(PyTuple_GetItem(t, 0)) ||
        !PyString_Check(PyTuple_GetItem(t, 1))) {
        PyErr_SetString(PyExc_ValueError, "tuples must contain strings");
        return NULL;
    }
    string ParameterName = PyString_AsString(s);
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
      cout << "Parameter type not recognized" << endl;
      exit(EXIT_FAILURE);
    }
  }
}

%typemap(freearg) (Teuchos::ParameterList& ParameterList)
{
   if ($1) delete($1);
}
// Amesos interface includes
%include "Amesos_config.h"
%include "Amesos.h"
%include "Amesos_BaseSolver.h"
#ifdef HAVE_AMESOS_LAPACK
%include "Amesos_Lapack.h"
#endif
#ifdef HAVE_AMESOS_KLU
%include "Amesos_Klu.h"
#endif
#ifdef HAVE_AMESOS_UMFPACK
%include "Amesos_Umfpack.h"
#endif
#ifdef HAVE_AMESOS_SUPERLU
%include "Amesos_Superlu.h"
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
%include "Amesos_Superludist.h"
#endif
#ifdef HAVE_AMESOS_MUMPS
%include "Amesos_Mumps.h"
#endif

// Extensions for Amesos
%extend Amesos_BaseSolver 
{
  string __str__() {
    stringstream os;
    os << "*** Amesos_BaseSolver ***";
    return os.str();
  }

  void __del__()
  {
    delete self;
  }
}


