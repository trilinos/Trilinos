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

// Do the local includes first to make sure that Python.h is included before
// any other standard headers, cf. <http://docs.python.org/2/c-api/intro.html#includes>.
#include "PyTrilinos_Util.h"
#include "PyTrilinos_PythonException.h"
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
