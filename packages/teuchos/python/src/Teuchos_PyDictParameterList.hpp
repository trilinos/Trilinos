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
#include "Teuchos_FILEstream.hpp"

namespace Teuchos {

  // Forward declarations
  bool setPythonParameter(ParameterList &, const std::string &, PyObject *);
  PyObject * getPythonParameter(const ParameterList &, const std::string &);
  bool isEquivalent(PyObject *, const ParameterList &);

  class PyDictParameterList : public ParameterList {

  public:

    // -------------------------------------------------------------------------
    // Default constructor returns PyDictParameterList with empty
    // ParameterList and empty python dictionary
    PyDictParameterList() :
      ParameterList()
    {
      // Create an empty python dictionary
      __pyDict__ = PyDict_New();
    }

    // -------------------------------------------------------------------------
    // String constructor returns PyDictParameterList with named empty
    // ParameterList and empty python dictionary
    PyDictParameterList(const string & name) :
      ParameterList(name)
    {
      // Create an empty python dictionary
      __pyDict__ = PyDict_New();
    }

    // -------------------------------------------------------------------------
    // Python dictionary constructor returns PyDictParameterList with
    // synchronized ParameterList and owned reference to provided
    // dictionary
    PyDictParameterList(PyObject * dict, string lName = string("ANONYMOUS")) :
      ParameterList(lName)
    {
      PyObject * key   = NULL;
      PyObject * value = NULL;
      int        pos   = 0;
      char     * name;

      // Check for a python dictionary
      if (!PyDict_Check(dict)) {
	PyErr_SetString(PyExc_TypeError, "Expecting a python dictionary");
	__pyDict__ = NULL;
	goto fail;
      }

      // Take ownership of a reference to the python dictionary
      Py_INCREF(dict);
      __pyDict__ = dict;

      // Iterate over the python dictionary entries and create
      // corresponding ParameterList entries
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
    // Copy constructor returns copy of source PyDictParameterList
    PyDictParameterList(const PyDictParameterList & source) :
      ParameterList(source)
    {
      // Copy the python dictionary
      __pyDict__ = PyDict_Copy(source.__pyDict__);
    }

    // -------------------------------------------------------------------------
    // ParameterList constructor returns PyDictParameterList with copy
    // of provided ParameterList and synchronized python dictionary
    PyDictParameterList(const ParameterList & source) :
      ParameterList(source)
    {
      // Create an empty python dictionary
      __pyDict__ = PyDict_New();

      // Iterate over the ParameterList entries and create
      // corresponding python dictionary entries
      for (ConstIterator i = begin(); i != end(); i++) {
	const string & name  = this->name(i);
	PyObject     * value = getPythonParameter(*this, name);

	// Unsupported python type
	if (value == NULL) {
	  PyErr_SetString(PyExc_TypeError,
			  "Source ParameterList has entry of unsupported python type");
	  goto fail;
	}
	// Parameter not found
	else if (value == Py_None) {
	  PyErr_Format(PyExc_KeyError, "'%s'", name.c_str());
	  goto fail;
	}
	// Parameter found and type supported
	if (PyDict_SetItemString(__pyDict__, name.c_str(), value) == -1)
	  goto fail;
      }
      return;
    fail:
      Py_DECREF(__pyDict__);
    }

    // -------------------------------------------------------------------------
    ~PyDictParameterList()
    {
      // Release ownership of the reference to the python dictionary.
      // If the constructor failed, the pointer could be NULL
      Py_XDECREF(__pyDict__);
    }

    // -------------------------------------------------------------------------
    // Assignment operator
    PyDictParameterList & operator=(const PyDictParameterList &source)
    {
      // Do nothing if self-assignment is requested
      if (&source == this) return *this;
      // Perform base class copy
      ParameterList::operator=(source);
      // Copy the python dictionary (do not take ownership of old reference)
      Py_XDECREF(__pyDict__);
      __pyDict__ = PyDict_Copy(source.__pyDict__);
      return *this;
    }

    /////////////////////////////////////////////////
    // Overloaded or overridden ParameterList methods
    /////////////////////////////////////////////////

    // -------------------------------------------------------------------------
    void set(const string &name, PyObject * value)
    {
      // Set the ParameterList parameter
      TEST_FOR_EXCEPTION(!setPythonParameter(*this,name,value), std::runtime_error, 
			 "ParameterList value type not supported");
      // Set the python dictionary parameter
      PyDict_SetItemString(__pyDict__, name.c_str(), value);
    }

    // -------------------------------------------------------------------------
    void set(const string &name, ParameterList plist)
    {
      sublist(name) = plist;
    }

    // -------------------------------------------------------------------------
    void set(const string &name, PyDictParameterList pdplist)
    {
      // Set the ParameterList sublist
      ParameterList::sublist(name) = ParameterList(pdplist);
      // Set the python dictionary parameter
      PyDict_SetItemString(__pyDict__, name.c_str(), pdplist.__pyDict__);
    }


    // -------------------------------------------------------------------------
    ParameterList & sublist(const string &name)
    {
      bool needToUpdate = not isParameter(name);
      ParameterList & result = ParameterList::sublist(name);
      if (needToUpdate) {
	PyDictParameterList pdplist(result);
	PyDict_SetItemString(__pyDict__, name.c_str(), pdplist.__pyDict__);
      }
      return result;
    }

//   ///////////////////////////////////////////////////
//   // Python methods that provide PyDict-like behavior
//   ///////////////////////////////////////////////////

    // -------------------------------------------------------------------------
    int __cmp__(PyObject * ob) const
    {
      return PyObject_Compare(__pyDict__, ob);
    }

    // -------------------------------------------------------------------------
    int __cmp__(const PyDictParameterList & plist) const
    {
      return PyObject_Compare(__pyDict__, plist.__pyDict__);
    }

    // -------------------------------------------------------------------------
    int __contains__(const char * key) const
    {
      // There is a simpler form in python 2.4, but this should be
      // backward-compatible
      PyObject * keys   = PyDict_Keys(__pyDict__);
      PyObject * keyStr = PyString_FromString(key);
      int result = PySequence_Contains(keys,keyStr);
      Py_DECREF(keys  );
      Py_DECREF(keyStr);
      return result;
    }

//     // -------------------------------------------------------------------------
//     void __delitem__(const char * key)
//     {
      
//     }

    // -------------------------------------------------------------------------
    PyObject * __eq__(PyObject * ob) const
    {
      return PyObject_RichCompare(__pyDict__, ob, Py_EQ);
    }

    // -------------------------------------------------------------------------
    PyObject * __eq__(const PyDictParameterList & plist) const
    {
      return PyObject_RichCompare(__pyDict__, plist.__pyDict__, Py_EQ);
    }

    // -------------------------------------------------------------------------
    PyObject * __ge__(PyObject * ob) const
    {
      return PyObject_RichCompare(__pyDict__, ob, Py_GE);
    }

    // -------------------------------------------------------------------------
    PyObject * __ge__(const PyDictParameterList & plist) const
    {
      return PyObject_RichCompare(__pyDict__, plist.__pyDict__, Py_GE);
    }

//   // -------------------------------------------------------------------------
//   PyObject * __getattribute__(char * name) const

    // -------------------------------------------------------------------------
    PyObject * __getitem__(const char * key) const
    {
      PyObject * result = PyDict_GetItemString(__pyDict__, key);
      if (result == NULL) PyErr_Format(PyExc_KeyError,"'%s'",key);
      return result;
    }

    // -------------------------------------------------------------------------
    PyObject * __gt__(PyObject * ob) const
    {
      return PyObject_RichCompare(__pyDict__, ob, Py_GT);
    }

    // -------------------------------------------------------------------------
    PyObject * __gt__(const PyDictParameterList & plist) const
    {
      return PyObject_RichCompare(__pyDict__, plist.__pyDict__, Py_GT);
    }

//   // -------------------------------------------------------------------------
//   PyObject * __iter__() const

    // -------------------------------------------------------------------------
    PyObject * __le__(PyObject * ob) const
    {
      return PyObject_RichCompare(__pyDict__, ob, Py_LE);
    }

    // -------------------------------------------------------------------------
    PyObject * __le__(const PyDictParameterList & plist) const
    {
      return PyObject_RichCompare(__pyDict__, plist.__pyDict__, Py_LE);
    }

    // -------------------------------------------------------------------------
    int __len__() const
    {
      return PyDict_Size(__pyDict__);
    }

    // -------------------------------------------------------------------------
    PyObject * __lt__(PyObject * ob) const
    {
      return PyObject_RichCompare(__pyDict__, ob, Py_LT);
    }

    // -------------------------------------------------------------------------
    PyObject * __lt__(const PyDictParameterList & plist) const
    {
      return PyObject_RichCompare(__pyDict__, plist.__pyDict__, Py_LT);
    }

    // -------------------------------------------------------------------------
    PyObject * __ne__(PyObject * ob) const
    {
      return PyObject_RichCompare(__pyDict__, ob, Py_NE);
    }

    // -------------------------------------------------------------------------
    PyObject * __ne__(const PyDictParameterList & plist) const
    {
      return PyObject_RichCompare(__pyDict__, plist.__pyDict__, Py_NE);
    }

    // -------------------------------------------------------------------------
    PyObject * __repr__() const
    {
      return PyObject_Repr(__pyDict__);
    }

//   // -------------------------------------------------------------------------
//   int        __setitem__(const char * key, PyObject * value)

    // -------------------------------------------------------------------------
    PyObject * __str__() const
    {
      return PyObject_Str(__pyDict__);
    }

//   // -------------------------------------------------------------------------
//   void       clear()

    // -------------------------------------------------------------------------
    PyDictParameterList copy() const
    {
      return PyDictParameterList(*this);
    }

    // -------------------------------------------------------------------------
    int has_key(const char * key) const
    {
      // There is a simpler form in python 2.4, but this should be
      // backward-compatible
      PyObject * keys   = PyDict_Keys(__pyDict__);
      PyObject * keyStr = PyString_FromString(key);
      int result = PySequence_Contains(keys,keyStr);
      Py_DECREF(keys  );
      Py_DECREF(keyStr);
      return result;
    }

    // -------------------------------------------------------------------------
    PyObject * items() const
    {
      return PyDict_Items(__pyDict__);
    }

//   // -------------------------------------------------------------------------
//   PyObject * iteritems() const

//   // -------------------------------------------------------------------------
//   PyObject * iterkeys() const

//   // -------------------------------------------------------------------------
//   PyObject * itervalues() const

    // -------------------------------------------------------------------------
    PyObject * keys() const {
      return PyDict_Keys(__pyDict__);
    }

//   // -------------------------------------------------------------------------
//   PyObject * pop(const char * key, PyObject * default = Py_None)

//   // -------------------------------------------------------------------------
//   PyObject * popitem()

//   // -------------------------------------------------------------------------
//   void       setdefault(const char * key, PyObject * value = Py_None)

//   // -------------------------------------------------------------------------
//   void       update(PyObject * dict)

//   // -------------------------------------------------------------------------
//   void       update(const PyDictParameterList & plist)

    // -------------------------------------------------------------------------
    PyObject * values() const
    {
      return PyDict_Values(__pyDict__);
    }

//   //////////////////////////////////////////////////
//   // Special methods specific to PyDictParameterList
//   //////////////////////////////////////////////////

    // -------------------------------------------------------------------------
    PyObject * dict() const
    {
      Py_INCREF(__pyDict__);
      return __pyDict__;
    }

    // -------------------------------------------------------------------------
    bool isSynchronized() const
    {
      return isEquivalent(__pyDict__, *this);
    }

  private:

    PyObject * __pyDict__;

  };   // class PyDictParameterList

}    // namespace Teuchos

#endif    // TEUCHOS_PYDICTPARAMETERLIST
