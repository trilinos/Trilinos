#ifndef PYTEUCHOS_UTILS_H
#define PYTEUCHOS_UTILS_H

#include "Teuchos_ParameterList.hpp"
#include "Python.h"

// Creates a newly allocated Teuchos parameter list from the input
// object, which must be a Python dictionary.
//
// "bool", "int", "double" and "string" are automatically recognized.
// Other types can be defined here as tuples. For example, the Python
// dictionary can be something like:
// List = {
//   "double parameter": 12.0,
//   "int parameter"   : 12,
//   "string parameter": "12"
// }
//
// \author Marzio Sala, SNL 9215
//
// \date Last modified on 08-Aug-05

Teuchos::ParameterList* CreateList(PyObject* obj)
{
  int i;
  Teuchos::ParameterList* List;
  if (!PyDict_Check(obj)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a dictionary");
    return NULL;
  }
  List = new Teuchos::ParameterList;

  int size = PyDict_Size(obj);
  PyObject* Keys = PyDict_Keys(obj);
  PyObject* Values = PyDict_Values(obj);

  for (i = 0; i < size ; i++) 
  {
    PyObject *s = PyList_GetItem(Keys,i);
    PyObject *t = PyList_GetItem(Values,i);

    // Get the parameter name
    if (!PyString_Check(s)) {
        PyErr_SetString(PyExc_ValueError, "Dictionary keys must be strings");
        return NULL;
    }
    string ParameterName = PyString_AsString(s);

    // now parse for the parameter value and type
    // This can be a "int", "double", "string", or a tuple
    // for more general types

    if (PyBool_Check(t)) 
    {
      if (t == Py_True)
        List->set(ParameterName, true);
      else
        List->set(ParameterName, false);
    }
    else if (PyInt_Check(t)) 
    {
      int ParameterValue = PyInt_AsLong(t);
      List->set(ParameterName, ParameterValue);
    }
    else if (PyFloat_Check(t)) 
    {
      double ParameterValue = PyFloat_AsDouble(t);
      List->set(ParameterName, ParameterValue);
    }
    else if (PyString_Check(t)) 
    {
      string ParameterValue = PyString_AsString(t);
      List->set(ParameterName, ParameterValue);
    }
    else if (PyTuple_Check(t)) 
    {
      if (!PyString_Check(PyTuple_GetItem(t, 0)) ||
          !PyString_Check(PyTuple_GetItem(t, 1))) {
        PyErr_SetString(PyExc_ValueError, "tuples must contain strings");
        return NULL;
      }
      string ParameterType = PyString_AsString(PyTuple_GetItem(t, 0));
      string ParameterValue = PyString_AsString(PyTuple_GetItem(t, 1));
      if (ParameterType == "bool") 
      {
        if (ParameterValue == "true")
          List->set(ParameterName, true);
        else
          List->set(ParameterName, false);
      }
      else if (ParameterType == "int") 
      {
        List->set(ParameterName, (int)atoi(ParameterValue.c_str()));
      }
      else if (ParameterType == "double") 
      {
        List->set(ParameterName, (double)atof(ParameterValue.c_str()));
      }
      else if (ParameterType == "string") 
      {
        List->set(ParameterName, string(ParameterValue));
      }
      else 
      {
        PyErr_SetString(PyExc_ValueError, "type in tuple not recognized");
        return NULL;
      }
    }
    else
    {
      PyErr_SetString(PyExc_ValueError, "Type in list not recognized");
      return NULL;
    }
  }

  return(List);
}
#endif
