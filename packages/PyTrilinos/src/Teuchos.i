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

// This documentation string will be the python help facility help
// string
%define DOCSTRING
"The Teuchos module allows access to The Trilinos package
Teuchos.  Use the python help() facility for local documentation
on classes and methods, or see the on-line documentation for more
in-depth information."
%enddef

// Define the module name, its package and documentation string
%module(package      = "PyTrilinos",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = DOCSTRING) Teuchos

// SWIG does not support wrapping nested classes.  We will %import the
// Teuchos::any class (ie, tell swig about it, but not wrap it), which
// has nested classes.  To suppress the swig warning that would
// otherwise result, we use the following:
#pragma SWIG nowarn=312

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
#include "Teuchos_PythonParameter.h"
%}

// Auto-documentation feature
%feature("autodoc", "1");

// C++ STL support.  If the wrapped class uses standard template
// library containers, the following %include wraps the containers
// and makes certain conversions seamless, such as between std::string
// and python strings.
%include "stl.i"
namespace std {
  class logic_error;
  class runtime_error;
}
using namespace std;

// Define macro for handling exceptions thrown by Teuchos methods and
// constructors
%define TEUCHOS_EXCEPTION(className,methodName)
%exception Teuchos::className::methodName {
  try {
    $action
    if (PyErr_Occurred()) SWIG_fail;
  }
  catch(Teuchos::Exceptions::InvalidParameterType e) {
    PyErr_SetString(PyExc_TypeError, e.what());
    SWIG_fail;
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

// General ignore directives
%ignore *::operator=;
%ignore *::print;

//Teuchos imports
%import "Teuchos_TypeNameTraits.hpp"
%import "Teuchos_RefCountPtrDecl.hpp"
%import "Teuchos_any.hpp"
%import "Teuchos_ParameterEntry.hpp"
%import "Teuchos_PythonParameter.h"

/////////////////////////////
// Teuchos_Version support //
/////////////////////////////
%include "Teuchos_Version.hpp"
%pythoncode %{
  __version__ = Teuchos_Version().split()[2]
%}

///////////////////////////////////
// Teuchos_ParameterList support //
///////////////////////////////////
TEUCHOS_EXCEPTION(ParameterList,ParameterList)
TEUCHOS_EXCEPTION(ParameterList,set)
TEUCHOS_EXCEPTION(ParameterList,setParameters)
TEUCHOS_EXCEPTION(ParameterList,get)
TEUCHOS_EXCEPTION(ParameterList,sublist)
TEUCHOS_EXCEPTION(ParameterList,type)
TEUCHOS_EXCEPTION(ParameterList,__setitem__)
TEUCHOS_EXCEPTION(ParameterList,update)
// There are a lot of extensions to the Teuchos::ParameterList class,
// so I put them all in their own file
%include "Teuchos_ParameterList_ext.i"
%ignore Teuchos::ParameterList::set;
%ignore Teuchos::ParameterList::setEntry;
%ignore Teuchos::ParameterList::get;
%ignore Teuchos::ParameterList::getPtr;
%ignore Teuchos::ParameterList::getEntryPtr;
%ignore Teuchos::ParameterList::sublist(const string &) const;
%ignore Teuchos::ParameterList::isType(const string &) const;
%ignore Teuchos::ParameterList::isType(const string &, any*) const;
%ignore Teuchos::ParameterList::unused(ostream &) const;
%ignore Teuchos::ParameterList::begin() const;
%ignore Teuchos::ParameterList::end() const;
%ignore Teuchos::ParameterList::entry(ConstIterator) const;
%ignore Teuchos::ParameterList::name(ConstIterator) const;
%include "Teuchos_ParameterList.hpp"

////////////////////////////////////////////////////////////////////////

// Typemaps.  These are generally intended for other packages that
// import this interface file.

// These typemaps allow C++ methods that expect ParameterList&
// arguments to have python wrappers that accept ParameterList or
// python dictionary arguments.

%typemap(in) Teuchos::ParameterList & (void * argp = 0, int res = 0, bool cleanup=false)
{
  if (PyDict_Check($input)) {
    $1 = Teuchos::pyDictToNewParameterList($input);
    if ($1 == NULL) SWIG_fail;
    cleanup = true;
  }

  else {
    res = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
    if (!SWIG_IsOK(res)) {
      %argument_fail(res, "$type", $symname, $argnum);
    }
    $1 = %reinterpret_cast(argp, $ltype);
  }
}

%typecheck(200) Teuchos::ParameterList &
{
  // Accept PyDicts or ParameterLists
  void * argp = NULL;
  $1 = PyDict_Check($input) ? 1 : 0;
  if (!$1) if (SWIG_CheckState(SWIG_Python_ConvertPtr($input, &argp,
						      $1_descriptor, 0))) $1 = 1;
}

%typemap(freearg) Teuchos::ParameterList &
{
  if (cleanup$argnum && $1) delete $1;
}

// Extend the %extend_smart_pointer macro to handle derived classes.
// This has been specialized for RefCountPtr, because the constructor
// for creating a RefCountPtr has an additional argument: a boolean
// false
%define %extend_RefCountPtr(Type, SmartPtrType...)
%typemap(in, noblock=1) const SWIGTYPE & SMARTPOINTER (void* argp = 0, int res = 0) {
  res = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor(Type*), %convertptr_flags);
    if (!SWIG_IsOK(res)) {
      %argument_fail(res, "$type", $symname, $argnum);
    }
    if (!argp) { %argument_nullref("$type", $symname, $argnum); }
    $1 = new $*ltype ( %reinterpret_cast( argp, Type* ), false );
  }
  else {
    if (!argp) { %argument_nullref("$type", $symname, $argnum); }
    $1 = %reinterpret_cast(argp, $ltype);
  }
}

%typecheck(1200) const SWIGTYPE & SMARTPOINTER {
  static void * argp = 0;
  $1 = SWIG_CheckState(SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags)) ? 1 : 0;
  if (!$1)
    $1 = SWIG_CheckState(SWIG_ConvertPtr($input, &argp, $descriptor(Type *),
					 %convertptr_flags)) ? 1 : 0;
}

%typemap(freearg) const SWIGTYPE & SMARTPOINTER {
  delete $1;
}

%extend_smart_pointer(SmartPtrType)

%enddef

// These typemap macros allow developers to generate typemaps for any
// classes that are wrapped in RefCountPtr as function or method
// arguments.
%define TEUCHOS_RCP_TYPEMAPS(Type)

%extend_RefCountPtr(Type, Teuchos::RefCountPtr< Type >)
%template()               Teuchos::RefCountPtr< Type >;

%extend_RefCountPtr(Type, Teuchos::RefCountPtr< const Type >)
%template()               Teuchos::RefCountPtr< const Type >;

%enddef

// The following directives are for the special case of a
// Teuchos::RefCountPtr that points to a Teuchos::ParameterList

// Apply the RefCountPtr typemap macros to ParameterLists
%ignore Teuchos::RefCountPtr< Teuchos::ParameterList >::get() const;
TEUCHOS_RCP_TYPEMAPS(Teuchos::ParameterList)
