// @HEADER
// ***********************************************************************
//
//           PyTrilinos.Teuchos: Python Interface to Teuchos
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

#include "PyDictParameterList.hpp"

// -------------------------------------------------------------------------
PyDictParameterList::PyDictParameterList() :
  Teuchos::ParameterList()
{
  __pyDict__ = PyDict_New();
}

// -------------------------------------------------------------------------
PyDictParameterList::PyDictParameterList(PyObject * dict) :
  Teuchos::ParameterList()
{
  PyObject *key    = NULL;
  PyObject *value  = NULL;
  int       result = 0;
  int       pos    = 0;
  if (!PyDict_Check(dict)) {
    PyErr_SetString(PyExc_TypeError, "Expecting type dict");
    __pyDict__ = NULL;
    goto fail;
  }
  Py_INCREF(dict);
  __pyDict__ = dict;
  while (PyDict_Next(__pyDict__, &pos, &key, &value)) {
    result = __setitem__(key,value);
    if (result < 0) goto fail;
  }
  return;
 fail:
  Py_XDECREF(__pyDict__);
}

// -------------------------------------------------------------------------
PyDictParameterList::PyDictParameterList(const PyDictParameterList & source) :
  Teuchos::ParameterList(source)
{
  __pyDict__ = PyDict_New();
}

// -------------------------------------------------------------------------
PyDictParameterList::~PyDictParameterList() {
  Py_DECREF(__pyDict__);
}


// -------------------------------------------------------------------------
int PyDictParameterList::__cmp__(PyObject * dict) const {

}

// -------------------------------------------------------------------------
int PyDictParameterList::__cmp__(const PyDictParameterList) const {

}

// -------------------------------------------------------------------------
bool PyDictParameterList::__contains__(const char * key) const {

}

// -------------------------------------------------------------------------
void PyDictParameterList::__delitem__(const char * key) {

}

// -------------------------------------------------------------------------
bool PyDictParameterList::__eq__(PyObject * dict) const {

}

// -------------------------------------------------------------------------
bool PyDictParameterList::__eq__(const PyDictParameterList) const {

}

// -------------------------------------------------------------------------
bool PyDictParameterList::__ge__(PyObject * dict) const {

}

// -------------------------------------------------------------------------
bool PyDictParameterList::__ge__(const PyDictParameterList) const {

}

// -------------------------------------------------------------------------
PyObject * PyDictParameterList::__getattribute__(char * name) const {

}

// -------------------------------------------------------------------------
PyObject * PyDictParameterList::__getitem__(const char * key) const {

}

// -------------------------------------------------------------------------
bool PyDictParameterList::__gt__(PyObject * dict) const {

}

// -------------------------------------------------------------------------
bool PyDictParameterList::__gt__(const PyDictParameterList) const {

}

// -------------------------------------------------------------------------
PyObject * PyDictParameterList::__iter__() const {

}

// -------------------------------------------------------------------------
bool PyDictParameterList::__le__(PyObject * dict) const {

}

// -------------------------------------------------------------------------
bool PyDictParameterList::__le__(const PyDictParameterList) const {

}

// -------------------------------------------------------------------------
int  PyDictParameterList::__len__() const {

}

// -------------------------------------------------------------------------
bool PyDictParameterList::__lt__(PyObject * dict) const {

}

// -------------------------------------------------------------------------
bool PyDictParameterList::__lt__(const PyDictParameterList) const {

}

// -------------------------------------------------------------------------
bool PyDictParameterList::__ne__(PyObject * dict) const {

}

// -------------------------------------------------------------------------
bool PyDictParameterList::__ne__(const PyDictParameterList) const {

}

// -------------------------------------------------------------------------
char * PyDictParameterList::__repr__() const {

}

// -------------------------------------------------------------------------
int PyDictParameterList::__setitem__(PyObject * key, PyObject * value) {

  string name;
  int    result;

  // Check that key is a string
  if (!PyString_Check(key)) {
    PyErr_SetString(PyExc_TypeError, "Expecting string keys in dict");
    goto fail;
  }
  name = PyString_AsString(key);

  // Boolean values
  if (PyBool_Check(value)) {
    if (value == Py_True) set(name,true );
    else                  set(name,false);
  }

  // Integer values
  else if (PyInt_Check(value)) {
    set(name, (int)PyInt_AsLong(value));
  }

  // Floating point values
  else if (PyFloat_Check(value)) {
    set(name, PyFloat_AsDouble(value));
  }

  // String values
  else if (PyString_Check(value)) {
    set(name, PyString_AsString(value));
  }

  // Dictionary values
  else if (PyDict_Check(value)) {
    PyDictParameterList * newList = new PyDictParameterList(value);
    set(name, newList);
  }

  // Unsupported value types
  else {
    PyErr_SetString(PyExc_TypeError, "Dictionary value type not supported");
    goto fail;
  }

  // Update the __pyDict__ object
  result = PyDict_SetItem(__pyDict, key, value);
  if (result < 0) goto fail;

  return 0;

 fail:
  return -1;
}

// -------------------------------------------------------------------------
char * PyDictParameterList::__str__() const {

}

// -------------------------------------------------------------------------
void PyDictParameterList::clear() {

}

// -------------------------------------------------------------------------
PyObject * PyDictParameterList::copy() const {

}

// -------------------------------------------------------------------------
PyObject * PyDictParameterList::get(const char * key, PyObject * value = Py_None) const {

}

// -------------------------------------------------------------------------
bool PyDictParameterList::has_key(const char * key) const {

}

// -------------------------------------------------------------------------
PyObject * PyDictParameterList::items() const {

}

// -------------------------------------------------------------------------
PyObject * PyDictParameterList::iteritems() const {

}

// -------------------------------------------------------------------------
PyObject * PyDictParameterList::iterkeys() const {

}

// -------------------------------------------------------------------------
PyObject * PyDictParameterList::itervalues() const {

}

// -------------------------------------------------------------------------
PyObject * PyDictParameterList::keys() const {

}

// -------------------------------------------------------------------------
PyObject * PyDictParameterList::pop(const char * key, PyObject * default = Py_None) {

}

// -------------------------------------------------------------------------
PyObject * PyDictParameterList::popitem() {

}

// -------------------------------------------------------------------------
void PyDictParameterList::setdefault(const char * key, PyObject * value = Py_None) {

}

// -------------------------------------------------------------------------
void PyDictParameterList::update(PyObject * dict) {

}

// -------------------------------------------------------------------------
void PyDictParameterList::update(const PyDictParameterList &) {

}

// -------------------------------------------------------------------------
PyObject * PyDictParameterList::values() const {

}
