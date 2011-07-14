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

#ifndef PYTRILINOS_UTIL_H
#define PYTRILINOS_UTIL_H

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

#endif // PYTRILINOS_UTIL_H
