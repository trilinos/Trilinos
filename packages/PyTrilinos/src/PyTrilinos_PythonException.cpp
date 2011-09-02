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

#include "PyTrilinos_PythonException.h"

namespace PyTrilinos
{

PythonException::PythonException() :
  std::runtime_error("")
{
  // If there is a python error, this python API function will place
  // the error information into the data members, with an owned
  // reference to each, and then clear the python error.  If there is
  // no python error, the data members will be set to NULL.
  PyErr_Fetch(&errorType, &errorValue, &errorTraceback);

  // If errorType is NULL, then the PythonException class is being
  // used improperly
  if (errorType == NULL)
  {
    Py_INCREF(PyExc_RuntimeError);
    errorType  = PyExc_RuntimeError;
    errorValue = PyString_FromString(
      "A C++ PythonException was thrown without a corresponding python error.\n"
      "Please contact a PyTrilinos developer with the details.");
  }

  // Obtain the python string representation of errorValue, convert it
  // to a char*, and decrement the reference to the python string
  // representation.
  PyObject * msg = PyObject_Str(errorValue);
  errorMsg       = PyString_AsString(msg);
  Py_DECREF(msg);
}

PythonException::~PythonException() throw()
{
  // Release all owned references
  Py_XDECREF(errorType     );
  Py_XDECREF(errorValue    );
  Py_XDECREF(errorTraceback);
}

const char * PythonException::what() const throw()
{
  return errorMsg;
}

void PythonException::restore()
{
  // This python API function resets the python error values and
  // reclaims the references to the python objects.  Therefore, set
  // the error data members to NULL.
  PyErr_Restore(errorType, errorValue, errorTraceback);
  errorType      = NULL;
  errorValue     = NULL;
  errorTraceback = NULL;
  errorMsg       = NULL;
}

}  // Namespace PyTrilinos
