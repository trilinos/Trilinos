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

// Do the local includes first to make sure that Python.h is included before
// any other standard headers, cf. <http://docs.python.org/2/c-api/intro.html#includes>.
#include "PyTrilinos_Util.hpp"
#include "PyTrilinos_PythonException.hpp"
#include "swigpyrun.h"

// System includes
#include <algorithm>

////////////////////////////////////////////////////////////////////////

namespace PyTrilinos
{

////////////////////////////////////////////////////////////////////////

PyObject * getObjectFromModule(char * modName, CONST char * objName)
{
  PyObject * fromList = NULL;
  PyObject * module   = NULL;
  PyObject * object   = NULL;
  fromList = Py_BuildValue("[s]", objName);
  if (!fromList) goto fail;
  module = PyImport_ImportModuleEx(modName, NULL, NULL, fromList);
  if (!module) goto fail;
  object = PyObject_GetAttrString(module, objName);
  if (!object) goto fail;
  Py_DECREF(fromList);
  Py_DECREF(module  );
  return object;

  fail:
  Py_XDECREF(fromList);
  Py_XDECREF(module);
  return NULL;
}

////////////////////////////////////////////////////////////////////////

PyObject * getClassFromModule(char * modName, CONST char * clsName)
{
  PyObject * cls = NULL;
  cls = getObjectFromModule(modName, clsName);
  if (!cls) goto fail;
  if (!PyType_Check(cls))
  {
    PyErr_Format(PyExc_TypeError, "Object '%s' is not a class type", clsName);
    goto fail;
  }
  return cls;

  fail:
  Py_XDECREF(cls);
  return NULL;
}

////////////////////////////////////////////////////////////////////////

bool objectAttrIsNone(PyObject * object, CONST char * name)
{
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!value) throw PythonException();
  bool result = (value == Py_None);
  Py_DECREF(value);
  return result;
}

////////////////////////////////////////////////////////////////////////

bool objectAttrIsTrue(PyObject * object, CONST char * name)
{
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!value) throw PythonException();
  bool result = (value == Py_True);
  Py_DECREF(value);
  return result;
}

////////////////////////////////////////////////////////////////////////

bool getBoolObjectAttr(PyObject * object, CONST char * name)
{
  bool result;
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!value) throw PythonException();
  if (!PyBool_Check(value))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type boolean", name);
    throw PythonException();
  }
  if (value == Py_True) result = true;
  else                  result = false;
  Py_DECREF(value);
  return result;
}

////////////////////////////////////////////////////////////////////////

int getIntObjectAttr(PyObject * object, CONST char * name)
{
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!value) throw PythonException();
  int result = (int) PyInt_AsLong(value);
  if (PyErr_Occurred()) throw PythonException();
  Py_DECREF(value);
  return result;
}

////////////////////////////////////////////////////////////////////////

double getFloatObjectAttr(PyObject * object, CONST char * name)
{
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!value) throw PythonException();
  double result = PyFloat_AsDouble(value);
  if (PyErr_Occurred()) throw PythonException();
  Py_DECREF(value);
  return result;
}

////////////////////////////////////////////////////////////////////////

PyObject * getTupleObjectAttr(PyObject * object, CONST char * name)
{
  PyObject * result = PyObject_GetAttrString(object, name);
  if (!result) throw PythonException();
  if (!PyTuple_Check(result))
  {
    PyErr_Format(PyExc_TypeError, "Attribute '%s' is not of type tuple", name);
    Py_DECREF(result);
    throw PythonException();
  }
  return result;
}

////////////////////////////////////////////////////////////////////////

CONST char* getStringObjectAttr(PyObject * object, CONST char * name)
{
  PyObject * value = PyObject_GetAttrString(object, name);
  if (!value) throw PythonException();
  CONST char * result = PyString_AsString(value);
  if (PyErr_Occurred()) throw PythonException();
  Py_DECREF(value);
  return result;
}

////////////////////////////////////////////////////////////////////////

CONST char * getStringItemObjectAttr(PyObject * object, CONST char * name, int i)
{
  PyObject * tuple = getTupleObjectAttr(object, name);
  PyObject * item  = PyTuple_GetItem(tuple, i);
  Py_DECREF(tuple);
  if (!item) throw PythonException();
  CONST char * result = PyString_AsString(item);
  Py_DECREF(item);
  if (PyErr_Occurred()) throw PythonException();
  return result;
}

////////////////////////////////////////////////////////////////////////

}
