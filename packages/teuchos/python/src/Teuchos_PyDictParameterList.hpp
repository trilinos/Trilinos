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
	PyErr_SetString(PyExc_TypeError, "Expecting type dict");
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
      // Take ownership of a reference to the source
      // PyDictParameterList's python dictionary
      Py_INCREF(source.__pyDict__);
      __pyDict__ = source.__pyDict__;
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
      TEST_FOR_EXCEPTION(!setPythonParameter(*this,name,value), std::runtime_error, 
			 "ParameterList value type not supported");
    }

    // -------------------------------------------------------------------------
    void set(const string &name, ParameterList value)
    {
      ParameterList::sublist(name) = PyDictParameterList(value);
    }

    // -------------------------------------------------------------------------
    void set(const string &name, PyDictParameterList value)
    {
      sublist(name) = value;
    }

    // -------------------------------------------------------------------------
    PyDictParameterList & sublist(const string &name)
    {
      PyDictParameterList * pdpl_ptr;
      // If parameter exists, throw an exception for non-lists, or
      // return the parameter if it is a list
      if (isParameter(name)) {
	TEST_FOR_EXCEPTION(!isSublist(name), std::runtime_error,
			   "Parameter " << name << "exists, but is not a list");
	pdpl_ptr = getPtr<PyDictParameterList>(name);
	return (*pdpl_ptr);
      }
      // If the parameter does not exist, create a new, empty sublist,
      // insert it in the list and return its reference
      const PyDictParameterList newSublist(this->name()+std::string("->")+name);
      setEntry(name,ParameterEntry(newSublist,false,true));
      pdpl_ptr = getPtr<PyDictParameterList>(name);
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

    // -------------------------------------------------------------------------
    ostream& print(ostream& os, int indent, bool showTypes, bool showFlags) const
    {
      if (begin() == end()) {
	for (int j = 0; j < indent; j ++)
	  os << ' ';
	os << "[empty list]" << endl;
      }
      else { 
	// Print parameters first
	for (ConstIterator i = begin(); i != end(); ++i) 
	  {
	    const ParameterEntry &entry_i = entry(i);
	    if (entry_i.isList())
	      continue;
	    for (int j = 0; j < indent; ++j)
	      os << ' ';
	    os << name(i);
	    if(showTypes)
	      os << " : " << entry_i.getAny(false).typeName();
	    os << " = "; entry_i.leftshift(os,showFlags); os << endl;
	  }
	// Print sublists second
	for (ConstIterator i = begin(); i != end(); ++i) 
	  {
	    const ParameterEntry &entry_i = entry(i);
	    if (!entry_i.isList())
	      continue;
	    for (int j = 0; j < indent; j ++)
	      os << ' ';
	    os << name(i) << " -> " << endl;
	    getValue<PyDictParameterList>(entry_i).print(os, indent+2, showTypes, showFlags);
	  }
      }
      return os;
    }

    // -------------------------------------------------------------------------
    PyObject * _print(PyObject * pf=NULL, int indent=0, bool showTypes=false,
		      bool showFlags=true) {
      PyObject * returnObject = pf;
      // No arguments
      if (pf==NULL) {
	print(std::cout,indent,showTypes,showFlags);
	returnObject = Py_None;
      }

      // Given non-file pf argument
      else {
	if (!PyFile_Check(pf)) {
	  PyErr_SetString(PyExc_IOError, "_print() method expects a file object");
	  goto fail;
	}

	// Given file pf argument
	else {
	  std::FILE *f = PyFile_AsFile(pf);
	  FILEstream buffer(f);
	  std::ostream os(&buffer);
	  print(os,indent,showTypes,showFlags);
	}
      }
      Py_INCREF(returnObject);
      return returnObject;
    fail:
      return NULL;
    }

//    // -------------------------------------------------------------------------
//    bool isSublist(const string &name) const

//   ///////////////////////////////////////////////////
//   // Python methods that provide PyDict-like behavior
//   ///////////////////////////////////////////////////

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

//   //////////////////////////////////////////////////
//   // Special methods specific to PyDictParameterList
//   //////////////////////////////////////////////////

    // -------------------------------------------------------------------------
    bool isSynchronized() const
    {
      PyObject * key   = NULL;
      PyObject * value = NULL;
      PyObject * param = NULL;
      int        pos   = 0;
      string     name;
      // Check that all entries in ParameterList are also in the
      // python dictionary
      for (ConstIterator i = begin(); i != end(); i++) {
	name  = this->name(i);
	param = getPythonParameter(  *this     ,name.c_str());
	value = PyDict_GetItemString(__pyDict__,name.c_str());
	if (param == NULL)                                   goto fail;
	if (value == NULL)                                   goto fail;
	if (PyObject_RichCompareBool(param,value,Py_EQ) < 1) goto fail;
	Py_DECREF(param);
      }
      // Check that all entries in the python dictionary are also in
      // the ParameterList
      while (PyDict_Next(__pyDict__, &pos, &key, &value)) {
	if (!PyString_Check(key)) goto fail;
	name  = string(PyString_AsString(key));
	param = getPythonParameter(*this, name);
	if (!isParameter(name))                              goto fail;
	if (param == NULL     )                              goto fail;
	if (PyObject_RichCompareBool(param,value,Py_EQ) < 1) goto fail;
	Py_DECREF(param);
      }
      // All checks passed
      return true;

    fail:
      Py_XDECREF(param);
      return false;
    }

  private:

    PyObject * __pyDict__;

  };   // class PyDictParameterList

}    // namespace Teuchos

#endif    // TEUCHOS_PYDICTPARAMETERLIST
