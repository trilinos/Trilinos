i// @HEADER
// ***********************************************************************
//
//                 PyTrilinos: Rapid Prototyping Package
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

#include "Callback.h"

#include <iostream>

using std::cout;
using std::endl;

Callback::Callback()
  : mp_callback(0)
{
  // Nothing here yet
}

Callback::~Callback()
{
  Py_XDECREF(mp_callback);   /* Dispose of callback */
  mp_callback = 0;
}

PyObject * Callback::setFunction( PyObject * p_args)
{
    PyObject *p_result = NULL;
    PyObject *p_temp   = NULL;
    PyObject *p_func   = NULL;


    assert(0 != p_args && "Null argument passed to setFunction()");
    
    if (PyArg_ParseTuple(p_args, "O", &p_temp)) {
      // Assume that if this is a tuple, the item in the tuple is a
      // PyObject that is a pointer to a Python function. This is the
      // case when this function is called from Python.
      p_func = p_temp;
    } else {
      // Otherwise we assume that this function is directly passed a
      // PyObject that is a pointer to a Python function.  This is the
      // case when this function is called from C++.
      p_func = p_args;
    }

    if (!PyCallable_Check(p_func)) {
      PyErr_SetString(PyExc_TypeError,
		      "Function parameter must be callable");
      cout << "PyObject passed to function is not callable" << endl ;
      return NULL;
    }
    Py_XINCREF(p_func);          /* Add a reference to new callback */
    Py_XDECREF(mp_callback);     /* Dispose of previous callback    */
    mp_callback = p_func;        /* Remember new callback           */
    /* Boilerplate to return "None" */
    Py_INCREF(Py_None);
    p_result = Py_None;

    assert(0 != mp_callback && "Pointer to callback not set");
    return p_result;
}
 
PyObject * Callback::getFunction()
{
  assert (0 != mp_callback && "Callback function not yet assigned");
  return mp_callback;
}

const PyObject * Callback::getFunction() const
{
  assert (0 != mp_callback && "Callback function not yet assigned");
  return mp_callback;
}
