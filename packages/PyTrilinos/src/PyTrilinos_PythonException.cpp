// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "PyTrilinos_PythonException.hpp"

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
