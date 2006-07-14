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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_PYDICTPARAMETERLIST_H
#define TEUCHOS_PYDICTPARAMETERLIST_H

#include "Python.h"

#include "Teuchos_ParameterList.hpp"

namespace Teuchos {

  // Forward declarations
  bool setPythonParameter(ParameterList &, const std::string &, PyObject *);
  PyObject * getPythonParameter(const ParameterList &, const std::string &);

  class PyDictParameterList : public ParameterList {

  public:

    // -------------------------------------------------------------------------
    PyDictParameterList() :
      Teuchos::ParameterList()
    {
      __pyDict__ = PyDict_New();
    }

    // -------------------------------------------------------------------------
    PyDictParameterList(const string & name) :
      Teuchos::ParameterList(name)
    {
      __pyDict__ = PyDict_New();
    }

    // -------------------------------------------------------------------------
    PyDictParameterList(PyObject * dict) :
      Teuchos::ParameterList()
    {
      PyObject 	  * key   = NULL;
      PyObject 	  * value = NULL;
      int      	    pos   = 0;
      char        * name;
      if (!PyDict_Check(dict)) {
	PyErr_SetString(PyExc_TypeError, "Expecting type dict");
	__pyDict__ = NULL;
	goto fail;
      }
      Py_INCREF(dict);
      __pyDict__ = dict;
      while (PyDict_Next(__pyDict__, &pos, &key, &value)) {
	if (!PyString_Check(key)) {
	  PyErr_SetString(PyExc_KeyError,"Dictionary keys must be strings");
	  goto fail;
	}
	name = PyString_AsString(key);
	if (!setPythonParameter(*this,std::string(name),value)) {
	  PyErr_Format(PyExc_TypeError,"Parameter '%s' is of unsupported type",name);
	  goto fail;
	}
      }
      return;
    fail:
      Py_XDECREF(__pyDict__);
    }

    // -------------------------------------------------------------------------
    PyDictParameterList(const ParameterList & source) :
      Teuchos::ParameterList(source)
    {
      __pyDict__ = PyDict_New();
    }

//    // -------------------------------------------------------------------------
//    PyDictParameterList(const PyDictParameterList & source) :
//     Teuchos::ParameterList(source)
//    {
//      __pyDict__ = PyDict_New();
//    }

    // -------------------------------------------------------------------------
    ~PyDictParameterList()
    {
      Py_DECREF(__pyDict__);
    }

    // Overloaded / overridden ParameterList methods

    // -------------------------------------------------------------------------
    void set(const string &name, PyObject * value)
    {
      TEST_FOR_EXCEPTION(!setPythonParameter(*this,name,value), std::runtime_error, 
			 "ParameterList value type not supported");
    }

    // -------------------------------------------------------------------------
    void set(const string &name, ParameterList value)
    {
      sublist(name) = PyDictParameterList(value);
    }

    // -------------------------------------------------------------------------
    void set(const string &name, PyDictParameterList value)
    {
      sublist(name) = value;
    }

    // -------------------------------------------------------------------------
    PyDictParameterList & sublist(const string &name)
    {
      // If parameter exists, throw an exception for non-lists, or
      // return the parameter if it is a list
      if (isParameter(name)) {
	TEST_FOR_EXCEPTION(!isSublist(name), std::runtime_error,
			   "Parameter " << name << "exists, but is not a list");
	PyDictParameterList * pdpl_ptr = getPtr<PyDictParameterList>(name);
	return (*pdpl_ptr);
      }
      // If the parameter does not exist, create a new, empty sublist,
      // insert it in the list and return its reference
      const PyDictParameterList newSublist(this->name()+std::string("->")+name);
      setEntry(name,ParameterEntry(newSublist,false,true));
      PyDictParameterList * pdpl_ptr = getPtr<PyDictParameterList>(name);
      return (*pdpl_ptr);
    }

    // -------------------------------------------------------------------------
    const PyDictParameterList & sublist(const string &name) const {
      // If it does not exist, throw an error
      TEST_FOR_EXCEPTION( !isParameter(name), std::runtime_error,
			  " Parameter " << name << " is not a valid list!" );
      // If it does exist and is a list, return the list value.
      TEST_FOR_EXCEPTION( !isSublist(name), std::runtime_error,
			  " Parameter " << name << " is not a list!" );
      const PyDictParameterList * pdpl_ptr = getPtr<const PyDictParameterList>(name);
      return (*pdpl_ptr);
    }

//    // -------------------------------------------------------------------------
//    bool isSublist(const string &name) const

//   // Python methods that provide PyDict-like behavior

//   // -------------------------------------------------------------------------
//   int        __cmp__(PyObject * dict) const

//   // -------------------------------------------------------------------------
//   int        __cmp__(const PyDictParameterList) const

//   // -------------------------------------------------------------------------
//   bool       __contains__(const char * key) const

//   // -------------------------------------------------------------------------
//   void       __delitem__(const char * key)

//   // -------------------------------------------------------------------------
//   bool       __eq__(PyObject * dict) const

//   // -------------------------------------------------------------------------
//   bool       __eq__(const PyDictParameterList) const

//   // -------------------------------------------------------------------------
//   bool       __ge__(PyObject * dict) const

//   // -------------------------------------------------------------------------
//   bool       __ge__(const PyDictParameterList) const

//   // -------------------------------------------------------------------------
//   PyObject * __getattribute__(char * name) const

//   // -------------------------------------------------------------------------
//   PyObject * __getitem__(const char * key) const

//   // -------------------------------------------------------------------------
//   bool       __gt__(PyObject * dict) const

//   // -------------------------------------------------------------------------
//   bool       __gt__(const PyDictParameterList) const

//   // -------------------------------------------------------------------------
//   PyObject * __iter__() const

//   // -------------------------------------------------------------------------
//   bool       __le__(PyObject * dict) const

//   // -------------------------------------------------------------------------
//   bool       __le__(const PyDictParameterList) const

//   // -------------------------------------------------------------------------
//   int        __len__() const

//   // -------------------------------------------------------------------------
//   bool       __lt__(PyObject * dict) const

//   // -------------------------------------------------------------------------
//   bool       __lt__(const PyDictParameterList) const

//   // -------------------------------------------------------------------------
//   bool       __ne__(PyObject * dict) const

//   // -------------------------------------------------------------------------
//   bool       __ne__(const PyDictParameterList) const

//   // -------------------------------------------------------------------------
//   char *     __repr__() const

//   // -------------------------------------------------------------------------
//   int        __setitem__(const char * key, PyObject * value)

//   // -------------------------------------------------------------------------
//   char *     __str__() const

//   // -------------------------------------------------------------------------
//   void       clear()

//   // -------------------------------------------------------------------------
//   PyObject * copy() const

//   // -------------------------------------------------------------------------
//   PyObject * get(const char * key, PyObject * value = Py_None) const

//   // -------------------------------------------------------------------------
//   bool       has_key(const char * key) const

//   // -------------------------------------------------------------------------
//   PyObject * items() const

//   // -------------------------------------------------------------------------
//   PyObject * iteritems() const

//   // -------------------------------------------------------------------------
//   PyObject * iterkeys() const

//   // -------------------------------------------------------------------------
//   PyObject * itervalues() const

//   // -------------------------------------------------------------------------
//   PyObject * keys() const

//   // -------------------------------------------------------------------------
//   PyObject * pop(const char * key, PyObject * default = Py_None)

//   // -------------------------------------------------------------------------
//   PyObject * popitem()

//   // -------------------------------------------------------------------------
//   void       setdefault(const char * key, PyObject * value = Py_None)

//   // -------------------------------------------------------------------------
//   void       update(PyObject * dict)

//   // -------------------------------------------------------------------------
//   void       update(const PyDictParameterList &)

//   // -------------------------------------------------------------------------
//   PyObject * values() const

  private:

    PyObject * __pyDict__;

  };   // class PyDictParameterList

}    // namespace Teuchos

#endif    // TEUCHOS_PYDICTPARAMETERLIST
