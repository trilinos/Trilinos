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

#ifndef PYTRILINOS_PYTHONEXCEPTION_HPP
#define PYTRILINOS_PYTHONEXCEPTION_HPP

// The purpose of the PythonException class is to convert a python
// error to a C++ exception.  This is useful, especially in C++
// constructors that access the python/C API, when a python error is
// raised but a C++ exception should be thrown in order to ensure
// proper cleanup of a failed operation.

// When a python error is detected, or if the code sets a python
// error, this should be followed with
//
//     throw PyTrilinos::PythonException();
//
// The PythonException will extract the information it needs directly
// from the python API.

// Typically, the wrapper code should catch this exception and convert
// it back to a python error, using the restore() method.  Using SWIG,
// this can be done with the %exception directive:

// %exception
// {
//    try
//    {
//      $action
//    }
//    catch (PyTrilinos::PythonException &e)
//    {
//       e.restore();
//       SWIG_fail;
//    }
// }

#include "Python3Compat.hpp"
#include <stdexcept>

namespace PyTrilinos
{

class PythonException : public std::runtime_error
{
public:
  PythonException();
  ~PythonException() throw();
  const char * what() const throw();
  void restore();
private:
  PyObject * errorType;
  PyObject * errorValue;
  PyObject * errorTraceback;
  char     * errorMsg;
};

}  // Namespace PyTrilinos

#endif

#if defined(PyTrilinos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The PyTrilinos package is deprecated"
#endif
#endif

