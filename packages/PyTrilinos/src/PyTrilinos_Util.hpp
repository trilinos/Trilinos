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

#ifndef PYTRILINOS_UTIL_HPP
#define PYTRILINOS_UTIL_HPP

// Include the Python prototypes
#include "Python.h"

// Python developers improved the const-correctness of char* variables
// in the C-API with the advent of version 2.5.
#if PY_VERSION_HEX >= 0x02050000
#define CONST const
#else
#define CONST
#endif

namespace PyTrilinos
{

// Given the name of a python object, extract it from the named python
// module.  If it does not exist in the given module, return NULL.
PyObject * getObjectFromModule(char * modName, CONST char * objName);

// Given the name of a python class, extract it from the named python
// module.  If it does not exist in the given module, return NULL.
PyObject * getClassFromModule(char * modName, CONST char * clsName);

// Given a python object and an attribute name, return 'true' if the
// attribute exists and is python 'None'.  If the attribute does not
// exist, throw a PythonException.
bool objectAttrIsNone(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return 'true' if the
// attribute exists and is python 'True'.  If the attribute does not
// exist, throw a PythonException.
bool objectAttrIsTrue(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return the boolean
// value of the attribute.  If the attribute does not exist or the
// attribute is not a boolean, throw a PythonException.
bool getBoolObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return the integer
// value of the attribute.  If the attribute does not exist or the
// attribute is not an integer, throw a PythonException.
int getIntObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return the floating
// point (double) value of the attribute.  If the attribute does not
// exist or the attribute is not a float, throw a PythonException.
double getFloatObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return the tuple value
// of the attribute.  If the attribute does not exist or the attribute
// is not a tuple, throw a PythonException.
PyObject * getTupleObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return the string
// value of the attribute.  If the attribute does not exist or the
// attribute is not a string, throw a PythonException.
CONST char * getStringObjectAttr(PyObject * object, CONST char * name);

// Given a python object and an attribute name, return the string
// value of the i-th item of the attribute.  If the attribute does not
// exist or the attribute is not a sequence of strings, throw a
// PythonException.
CONST char * getStringItemObjectAttr(PyObject * object, CONST char * name, int i);

}  // Namespace PyTrilinos

#endif // PYTRILINOS_UTIL_HPP

#if defined(PyTrilinos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The PyTrilinos package is deprecated"
#endif
#endif

