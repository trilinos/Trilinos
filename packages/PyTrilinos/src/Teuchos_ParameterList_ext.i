// -*- c++ -*-

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

// The python implementation of Teuchos::ParameterList is augmented to
// "play nice" with python dictionaries, and to some extent, behave
// like a dictionary.  A constructor is added that takes a dictionary.
// In fact, any method that takes a ParameterList is extended to also
// take a dictionary.  Most dictionary methods and operators are added
// to the python implementation of ParameterList, so that it can be
// treated as a dictionary.  The exception to this is deletion
// methods, such as pop(), popitem() or __delitem__().  Since
// ParameterList does not support deletion, these methods are not
// implemented.

namespace Teuchos {

  %extend ParameterList {

    /******************************************************************/
    // Dictionary constructor
    ParameterList(PyObject * dict, string name = string("ANONYMOUS")) {
      Teuchos::ParameterList * plist = Teuchos::pyDictToNewParameterList(dict,
									 Teuchos::raiseError);
      if (plist == NULL) goto fail;

      plist->setName(name);
      return plist;
    fail:
      return NULL;
    }

    /******************************************************************/
    // Set method: accept only python objects as values
    PyObject * set(const string &name, PyObject *value) {
      if (!setPythonParameter(*self,name,value)) {
	PyErr_SetString(PyExc_TypeError, "ParameterList value type not supported");
	goto fail;
      }
      return Py_BuildValue("");
    fail:
      return NULL;
    }

    /******************************************************************/
    // SetParameters method, overloaded to accept a python dictionary
    ParameterList & setParameters(PyObject * dict) {
      if (!updateParameterListWithPyDict(dict,*self)) {
	PyErr_SetString(PyExc_ValueError, "ParameterList has values of unsupported type");
      }
      return *self;
    }

    /******************************************************************/
    // Get method: return entries as python objects
    PyObject * get(const string &name, PyObject * default_value=NULL) const {
      PyObject * value = getPythonParameter(*self, name);
      // Type not supported
      if (value == NULL) {
	PyErr_SetString(PyExc_TypeError, "ParameterList value type not supported");
	goto fail;
      }
      // Name not found
      else if (value == Py_None) {
	if (default_value == NULL) {
	  PyErr_Format(PyExc_KeyError, "'%s'", name.c_str());
	  goto fail;
	}
	Py_DECREF(value);
	Py_INCREF(default_value);
	return default_value;
      }
      // Type supported and name found
      else return value;
    fail:
      Py_XDECREF(value);
      return NULL;
    }

    /******************************************************************/
    // Print method: change name from "print" to "_print" because
    // "print" is a python keyword
    PyObject * _print(PyObject * pf=NULL, int indent=0, bool showTypes=false,
		      bool showFlags=true) {
      PyObject * returnObject = pf;
      // No arguments
      if (pf==NULL) {
	self->print(std::cout,indent,showTypes,showFlags);
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
	  Teuchos::FILEstream buffer(f);
	  std::ostream os(&buffer);
	  self->print(os,indent,showTypes,showFlags);
	}
      }
      Py_INCREF(returnObject);
      return returnObject;
    fail:
      return NULL;
    }

    /******************************************************************/
    // Unused method: take python file objects rather than ostreams
    PyObject * unused(PyObject * pf=NULL) {
      // No arguments
      if (pf==NULL) {
	self->unused(std::cout);
      }
      // Given non-file pf argument
      else {
	if (!PyFile_Check(pf)) {
	  PyErr_SetString(PyExc_IOError, "unused() method expects a file object");
	  goto fail;
	}
	// Given file pf argument
	else {
	  std::FILE *f = PyFile_AsFile(pf);
	  Teuchos::FILEstream buffer(f);
	  std::ostream os(&buffer);
	  self->unused(os);
	}
      }
      return Py_BuildValue("");
    fail:
      return NULL;
    }

    /******************************************************************/
    // Type method: return python type of requested parameter.  Replaces
    // templated C++ isType() methods.
    PyObject * type(const string & name) {
      PyObject * value = getPythonParameter(*self,name);
      // Type not supported
      if (value == NULL) {
	PyErr_SetString(PyExc_TypeError, "ParameterList value type not supported");
	goto fail;
      }
      // Name not found
      else if (value == Py_None) {
	PyErr_Format(PyExc_KeyError, "'%s'", name.c_str());
	goto fail;
      }
      // Name found and type supported
      return PyObject_Type(value);
    fail:
      return NULL;
    }

    //////////////////////////////////////////////
    // ** The following methods are added to ** //
    // ** give ParameterList a PyDict "feel" ** //
    //////////////////////////////////////////////

    /******************************************************************/
    // Comparison operators
    int __cmp__(PyObject * obj) const {
      PyObject * dict = Teuchos::parameterListToNewPyDict(*self,Teuchos::ignore);
      int result = 0;
      if (dict == NULL) goto fail;
      result = PyObject_Compare(dict,obj);
      Py_DECREF(dict);
      return result;
    fail:
      Py_XDECREF(dict);
      return -2;
    }

    int __cmp__(const ParameterList & plist) const {
      PyObject * dict1 = Teuchos::parameterListToNewPyDict(*self,Teuchos::ignore);
      PyObject * dict2 = Teuchos::parameterListToNewPyDict(plist,Teuchos::ignore);
      int result = 0;
      if (dict1 == NULL) goto fail;
      if (dict2 == NULL) goto fail;
      result = PyObject_Compare(dict1,dict2);
      Py_DECREF(dict1);
      Py_DECREF(dict2);
      return result;
    fail:
      Py_XDECREF(dict1);
      Py_XDECREF(dict2);
      return -2;
    }
    
    /******************************************************************/
    // Contains operator
    int __contains__(const string & name) const {
      PyObject * dict   = Teuchos::parameterListToNewPyDict(*self,Teuchos::ignore);
      PyObject * keys   = 0;
      PyObject * keyStr = 0;
      int result        = 0;
      if (dict == NULL) goto fail;
      // There is a simpler form in python 2.4, but this should be
      // backward-compatible
      keys   = PyDict_Keys(dict);
      keyStr = PyString_FromString(name.c_str());
      result = PySequence_Contains(keys,keyStr);
      Py_DECREF(dict  );
      Py_DECREF(keys  );
      Py_DECREF(keyStr);
      return result;
    fail:
      Py_XDECREF(dict);
      return result;
    }

    /******************************************************************/
    // Equals operators
    PyObject * __eq__(PyObject * obj) const {
      PyObject * dict   = Teuchos::parameterListToNewPyDict(*self,Teuchos::ignore);
      PyObject * result = 0;
      if (dict == NULL) goto fail;
      result = PyObject_RichCompare(dict,obj,Py_EQ);
      Py_DECREF(dict);
      return result;
    fail:
      return NULL;
    }

    PyObject * __eq__(const ParameterList & plist) const {
      PyObject * dict1  = Teuchos::parameterListToNewPyDict(*self,Teuchos::ignore);
      PyObject * dict2  = Teuchos::parameterListToNewPyDict(plist,Teuchos::ignore);
      PyObject * result = 0;
      if (dict1 == NULL) goto fail;
      if (dict2 == NULL) goto fail;
      result = PyObject_RichCompare(dict1,dict2,Py_EQ);
      Py_DECREF(dict1);
      Py_DECREF(dict2);
      return result;
    fail:
      Py_XDECREF(dict1);
      Py_XDECREF(dict2);
      return NULL;
    }

    /******************************************************************/
    // GetItem operator
    PyObject * __getitem__(const string & name) const {
      // I'm using SWIG's mangling scheme here
      // return Teuchos_ParameterList_get__SWIG_0(self,name);
      return Teuchos_ParameterList_get(self,name);
    }

    /******************************************************************/
    // __iter__ method
    PyObject * __iter__() const {
      PyObject * dict = Teuchos::parameterListToNewPyDict(*self,Teuchos::ignore);
      PyObject * iter = 0;
      if (dict == NULL) goto fail;
      iter = PyObject_GetIter(PyDict_Keys(dict));
      Py_DECREF(dict);
      return iter;
    fail:
      return NULL;
    }

    /******************************************************************/
    // Length operator
    int __len__() const {
      PyObject * dict = Teuchos::parameterListToNewPyDict(*self,Teuchos::ignore);
      int len = 0;
      if (dict == NULL) goto fail;
      len = PyDict_Size(dict);
      Py_DECREF(dict);
      return len;
    fail:
      return -1;
    }

    /******************************************************************/
    // Not equals operators
    PyObject * __ne__(PyObject * obj) const {
      PyObject * dict   = Teuchos::parameterListToNewPyDict(*self,Teuchos::ignore);
      PyObject * result = 0;
      if (dict == NULL) goto fail;
      result = PyObject_RichCompare(dict,obj,Py_NE);
      Py_DECREF(dict);
      return result;
    fail:
      return NULL;
    }

    PyObject * __ne__(const ParameterList & plist) const {
      PyObject * dict1  = Teuchos::parameterListToNewPyDict(*self,Teuchos::ignore);
      PyObject * dict2  = Teuchos::parameterListToNewPyDict(plist,Teuchos::ignore);
      PyObject * result = 0;
      if (dict1 == NULL) goto fail;
      if (dict2 == NULL) goto fail;
      result = PyObject_RichCompare(dict1,dict2,Py_NE);
      Py_DECREF(dict1);
      Py_DECREF(dict2);
      return result;
    fail:
      Py_XDECREF(dict1);
      Py_XDECREF(dict2);
      return NULL;
    }

    /******************************************************************/
    // SetItem operator
    void __setitem__(const string & name, PyObject * value) {
      // I'm using SWIG's mangling scheme here
      Teuchos_ParameterList_set(self,name,value);
    }

    /******************************************************************/
    // String representation method
    PyObject * __repr__() const {
      string reprStr;
      PyObject * dict    = Teuchos::parameterListToNewPyDict(*self,Teuchos::ignore);
      PyObject * dictStr = 0;
      PyObject * result = 0;
      if (dict == NULL) goto fail;
      dictStr = PyObject_Str(dict);
      reprStr = string("ParameterList(") + 
	        string(PyString_AsString(dictStr)) +
	        string(")");
      result = PyString_FromString(reprStr.c_str());
      Py_DECREF(dict   );
      Py_DECREF(dictStr);
      return result;
    fail:
      Py_XDECREF(dict);
      return NULL;
    }

    /******************************************************************/
    // String conversion method
    PyObject * __str__() const {
      PyObject * dict = Teuchos::parameterListToNewPyDict(*self,Teuchos::ignore);
      PyObject * str  = 0;
      if (dict == NULL) goto fail;
      str = PyObject_Str(dict);
      Py_DECREF(dict);
      return str;
    fail:
      return NULL;
    }

    /******************************************************************/
    // Has_key method
    int has_key(const string & name) const {
      PyObject * dict   = Teuchos::parameterListToNewPyDict(*self,Teuchos::ignore);
      PyObject * keys   = 0;
      PyObject * keyStr = 0;
      int result        = 0;
      if (dict == NULL) goto fail;
      // There is a simpler form in python 2.4, but this should be
      // backward-compatible
      keys   = PyDict_Keys(dict);
      keyStr = PyString_FromString(name.c_str());
      result = PySequence_Contains(keys,keyStr);
      Py_DECREF(dict  );
      Py_DECREF(keys  );
      Py_DECREF(keyStr);
      return result;
    fail:
      return -2;
    }

    /******************************************************************/
    // Items method
    PyObject * items() const {
      PyObject * dict   = Teuchos::parameterListToNewPyDict(*self,Teuchos::ignore);
      PyObject * result = 0;
      if (dict == NULL) goto fail;
      result = PyDict_Items(dict);
      Py_DECREF(dict);
      return result;
    fail:
      return NULL;
    }

    /******************************************************************/
    // Iteritems method
    PyObject * iteritems() const {
      PyObject * dict   = Teuchos::parameterListToNewPyDict(*self,Teuchos::ignore);
      PyObject * result = 0;
      if (dict == NULL) goto fail;
      result = PyObject_GetIter(PyDict_Items(dict));
      Py_DECREF(dict);
      return result;
    fail:
      return NULL;
    }

    /******************************************************************/
    // Iterkeys method
    PyObject * iterkeys() const {
      PyObject * dict   = Teuchos::parameterListToNewPyDict(*self,Teuchos::ignore);
      PyObject * result = 0;
      if (dict == NULL) goto fail;
      result = PyObject_GetIter(PyDict_Keys(dict));
      Py_DECREF(dict);
      return result;
    fail:
      return NULL;
    }

    /******************************************************************/
    // Itervalues method
    PyObject * itervalues() const {
      PyObject * dict   = Teuchos::parameterListToNewPyDict(*self,Teuchos::ignore);
      PyObject * result = 0;
      if (dict == NULL) goto fail;
      result = PyObject_GetIter(PyDict_Values(dict));
      Py_DECREF(dict);
      return result;
    fail:
      return NULL;
    }

    /******************************************************************/
    // Keys method
    PyObject * keys() const {
      PyObject * dict   = Teuchos::parameterListToNewPyDict(*self,Teuchos::ignore);
      PyObject * result = 0;
      if (dict == NULL) goto fail;
      result = PyDict_Keys(dict);
      Py_DECREF(dict);
      return result;
    fail:
      return NULL;
    }

    /******************************************************************/
    // Update methods
    void update(PyObject * dict, bool strict=true) {
      Teuchos::ResponseToIllegalParameters flag;
      if (strict) flag = Teuchos::raiseError;
      else        flag = Teuchos::storeNames;
      updateParameterListWithPyDict(dict,*self,flag);
    }

    void update(const ParameterList & plist) {
      self->setParameters(plist);
    }

    /******************************************************************/
    // Values method
    PyObject * values() const {
      PyObject * dict   = Teuchos::parameterListToNewPyDict(*self,Teuchos::ignore);
      PyObject * result = 0;
      if (dict == NULL) goto fail;
      result = PyDict_Values(dict);
      Py_DECREF(dict);
      return result;
    fail:
      return NULL;
    }

    /******************************************************************/
    // AsDict method: return a dictionary equivalent to the
    // ParameterList
    PyObject * asDict() const {
      return Teuchos::parameterListToNewPyDict(*self,Teuchos::storeNames);
    }

  }    // %extend ParameterList
}    // namespace Teuchos
