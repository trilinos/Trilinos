// -*- c++ -*-

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

// This documentation string will be the python help facility help
// string
%define DOCSTRING
"The Teuchos module allows access to The Trilinos package
Teuchos.  Use the python help() facility for local documentation
on classes and methods, or see the on-line documentation for more
in-depth information."
%enddef

// Define the module name, its package and documentation string
%module(package="PyTrilinos", docstring=DOCSTRING) Teuchos

// Code within the percent-bracket delimiters is copied verbatim to
// the C++ wrapper source file.  Anything that is %include-ed later
// needs to be #include-ed here.
%{
// System includes
#include <sstream>

// Teuchos includes
#include "Teuchos_FILEstream.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_RefCountPtrDecl.hpp"
#include "Teuchos_any.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_ParameterList.hpp"

// Teuchos python interface includes
#include "Teuchos_PythonParameter.hpp"
#include "Teuchos_PyDictParameterList.hpp"
%}

// Ignore directives.  Here we use them to prevent wrapping the Print
// methods.  Instead, we will later define __str__() methods for the
// classes, which is standard python technique for output.  The ignore
// directive is also good for eliminating warnings about methods
// python cannot wrap, such as operator=.
%ignore *::operator=;
%ignore *::print;
%ignore Teuchos::ParameterList::set(const string &, ParameterList);
%ignore Teuchos::ParameterList::set(const string &, const char[ ]);
%ignore Teuchos::ParameterList::setEntry(const string &, const ParameterEntry &);
%ignore Teuchos::ParameterList::get(const string &, const char[ ]);
%ignore Teuchos::ParameterList::getPtr(const string &);
%ignore Teuchos::ParameterList::getPtr(const string &) const;
%ignore Teuchos::ParameterList::getEntryPtr(const string &);
%ignore Teuchos::ParameterList::getEntryPtr(const string &) const;
%ignore Teuchos::ParameterList::sublist(const string &) const;
%ignore Teuchos::ParameterList::isType(const string &) const;
%ignore Teuchos::ParameterList::isType(const string &, any*) const;
%ignore Teuchos::ParameterList::unused(ostream &) const;
%ignore Teuchos::ParameterList::begin() const;
%ignore Teuchos::ParameterList::end() const;
%ignore Teuchos::ParameterList::entry(ConstIterator) const;
%ignore Teuchos::ParameterList::name(ConstIterator) const;

// Define macro for handling exceptions thrown by Teuchos methods and
// constructors
%define TEUCHOS_EXCEPTION(className,methodName)
%exception Teuchos::className::methodName {
  try {
    $action
    if (PyErr_Occurred()) SWIG_fail;
  }
  catch(Teuchos::Exceptions::InvalidParameter e) {
    PyErr_SetString(PyExc_KeyError, e.what());
    SWIG_fail;
  }
  catch(std::runtime_error e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    SWIG_fail;
  }
  catch(std::logic_error e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    SWIG_fail;
  }
}
%enddef

TEUCHOS_EXCEPTION(ParameterList,ParameterList)
TEUCHOS_EXCEPTION(ParameterList,set)
TEUCHOS_EXCEPTION(ParameterList,get)
TEUCHOS_EXCEPTION(ParameterList,sublist)
TEUCHOS_EXCEPTION(ParameterList,type)
TEUCHOS_EXCEPTION(PyDictParameterList,PyDictParameterList)
TEUCHOS_EXCEPTION(PyDictParameterList,set)
TEUCHOS_EXCEPTION(PyDictParameterList,get)
TEUCHOS_EXCEPTION(PyDictParameterList,setParameters)
TEUCHOS_EXCEPTION(PyDictParameterList,sublist)
TEUCHOS_EXCEPTION(PyDictParameterList,type)
TEUCHOS_EXCEPTION(PyDictParameterList,__contains__)
TEUCHOS_EXCEPTION(PyDictParameterList,__setitem__)
TEUCHOS_EXCEPTION(PyDictParameterList,has_key)
TEUCHOS_EXCEPTION(PyDictParameterList,update)

// Auto-documentation feature.  This ensures that calling help() on a
// wrapped method returns an argument list (or lists, in the case of
// overloaded methods).  The "1" option includes type information in
// the list.  While not as extensive as say, doxygen documentation,
// this is often enough to greatly increase the wrapper's usability.
%feature("autodoc", "1");

// C++ STL support.  If the wrapped class uses standard template
// library containers, the following %include-s wraps the containers
// and makes certain conversions seamless, such as between std::string
// and python strings.
%include "std_string.i"

// Teuchos interface includes
using namespace std;
%include "Teuchos_Version.hpp"
%import  "Teuchos_RefCountPtrDecl.hpp"
%import  "Teuchos_any.hpp"
%import  "Teuchos_ParameterEntry.hpp"
%import  "Teuchos_PythonParameter.hpp"
%include "Teuchos_ParameterList.hpp"
%include "Teuchos_PyDictParameterList.hpp"

// Extensions
%extend Teuchos::ParameterList {

  PyObject * set(const string &name, PyObject *value) {

    if (!setPythonParameter(*self,name,value)) {
      PyErr_SetString(PyExc_TypeError, "ParameterList value type not supported");
      goto fail;
    }

    return Py_BuildValue("");

  fail:
    return NULL;
  }

  PyObject * get(const string &name, PyObject * default_value=NULL) {

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

}

// Python code.  This code is added to the end of the python proxy
// file created by swig.  Here we set the namespace attribute
// "__version__" equal to the value returned by the
// Teuchos_Version function.
%pythoncode %{

  __version__ = Teuchos_Version().split()[2]

%}
