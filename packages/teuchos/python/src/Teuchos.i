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
#include "Teuchos_Version.hpp"
#include "Teuchos_RefCountPtrDecl.hpp"
#include "Teuchos_any.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_ParameterList.hpp"

// Teuchos::ParameterList can support parameters of any type, but the
// python wrappers need a subset of types supported a-priori.  Since
// this subset shows up in many places, the setPythonParameter()
// function is provided here that handles conversion from PyObjects to
// the supported types, assigning it to the specified ParameterList.
// It returns true if the conversion was made, false if the requested
// conversion is not supported.
  bool setPythonParameter(Teuchos::ParameterList & plist,
			  const std::string      & name,
			  PyObject               * value) {

    // Boolean values
    if (PyBool_Check(value)) {
      if (value == Py_True) plist.set(name,true );
      else                  plist.set(name,false);
    }

    // Integer values
    else if (PyInt_Check(value)) {
      plist.set(name, (int)PyInt_AsLong(value));
    }

    // Floating point values
    else if (PyFloat_Check(value)) {
      plist.set(name, PyFloat_AsDouble(value));
    }

    // String values
    else if (PyString_Check(value)) {
      plist.set(name, PyString_AsString(value));
    }

    // ParameterList values
    else {
      static swig_type_info * ty = SWIG_TypeQuery("Teuchos::ParameterList *");
      void * argp;
      int res = SWIG_ConvertPtr(value, &argp, ty, 0);
      if (SWIG_CheckState(res)) {
	Teuchos::ParameterList *arg = reinterpret_cast< Teuchos::ParameterList * >(argp);
	plist.set(name, &arg);
      }

      // Unsupported value types
      else {
	return false;
      }
    }
    return true;
  }

// The getPythonParameter() function is the get counterpart to
// setPythonParameter().  If the requested parameter name does not
// exist, None is returned (a type that is guaranteed not to be
// supported).  If the name exists and its type is supported, it is
// returned as a python object.  If the name exists, but the type is
// not supported, NULL is returned, to indicate an error.
  PyObject * getPythonParameter(const Teuchos::ParameterList & plist,
				const std::string            & name) {
    static swig_type_info * swig_TPL_ptr = SWIG_TypeQuery("Teuchos::ParameterList *");

    // Check for parameter existence
    if (!plist.isParameter(name)) return Py_BuildValue("");

    // Boolean parameter values
    if (Teuchos::isParameterType<bool>(plist,name)) {
      bool value = Teuchos::getParameter<bool>(plist,name);
      return PyBool_FromLong((long)value);
    }

    // Integer parameter values
    else if (Teuchos::isParameterType<int>(plist,name)) {
      int value = Teuchos::getParameter<int>(plist,name);
      return PyInt_FromLong((long)value);
    }

    // Double parameter values
    else if (Teuchos::isParameterType<double>(plist,name)) {
      double value = Teuchos::getParameter<double>(plist,name);
      return PyFloat_FromDouble(value);
    }

    // String parameter values
    else if (Teuchos::isParameterType<std::string>(plist,name)) {
      std::string value = Teuchos::getParameter<std::string>(plist,name);
      return PyString_FromString(value.c_str());
    }

    // Char* parameter values
    else if (Teuchos::isParameterType<char*>(plist,name)) {
      char* value = Teuchos::getParameter<char*>(plist,name);
      return PyString_FromString(value);
    }

    // ParameterList values
    else if (Teuchos::isParameterType<Teuchos::ParameterList>(plist,name)) {
      Teuchos::ParameterList value = Teuchos::getParameter<Teuchos::ParameterList>(plist,name);
      return SWIG_NewPointerObj((void*) &value, swig_TPL_ptr, 0);
    }

    // Unsupported type
    return NULL;
  }
%}

// Ignore directives.  Here we use them to prevent wrapping the Print
// methods.  Instead, we will later define __str__() methods for the
// classes, which is standard python technique for output.  The ignore
// directive is also good for eliminating warnings about methods
// python cannot wrap, such as operator=.
%ignore *::operator=;
%ignore *::print;
%ignore Teuchos::ParameterList::set(const string &, Teuchos::any);
%ignore Teuchos::ParameterList::set(const string &, Teuchos::ParameterList);
%ignore Teuchos::ParameterList::set(const string &, const char[ ]);
%ignore Teuchos::ParameterList::get(const string &, Teuchos::any);
%ignore Teuchos::ParameterList::get(const string &, const char[ ]);
%ignore Teuchos::ParameterList::getPtr(const string &);
%ignore Teuchos::ParameterList::getPtr(const string &) const;
%ignore Teuchos::ParameterList::getEntryPtr(const string &);
%ignore Teuchos::ParameterList::getEntryPtr(const string &) const;

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

// Teuchos interface includes.  Create a %include line for every
// header file with a prototype you want wrapped.  In this example,
// Teuchos_Version contains a function, Newp_Hello contains a
// class, and Newp_Jambo contains an optional class.
using namespace std;
%include "Teuchos_Version.hpp"
%include "Teuchos_RefCountPtrDecl.hpp"
%import  "Teuchos_any.hpp"
%include "Teuchos_ParameterEntry.hpp"
%include "Teuchos_ParameterList.hpp"

// Extensions.  The %extend directive allows you to extend classes to
// include additional methods.  In this example, we are adding a
// __str__() method, which is a standard python method for classes.
// When the python "print" command encounters a non-string object, it
// calls the str() function, which in turn looks for a __str__()
// method in the class.  Thus if hello is a Newp_Hello object, then
// "print hello" will return a string obtained from the C++ Print()
// method.
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
      Py_INCREF(default_value);
      return default_value;
    }

    // Type supported and name found
    else return value;

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
