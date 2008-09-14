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

%define %teuchos_docstring
"
PyTrilinos.Teuchos is the python interface to the Trilinos tools and
utilities package Teuchos:

    http://trilinos.sandia.gov/packages/teuchos

The purpose of Teuchos is to provide a number of utilities often
needed by numerical applications, but that are not necessarily
numerical by nature.  The python version of the Teuchos package
supports the following classes:

    * SerialComm              - Serial communicator
    * MpiComm                 - MPI communicator
    * DefaultComm             - Default communicator facility
    * ParameterList           - List of arbitrarily-typed values,
                                keyed by strings
    * XMLObject               - Object-oriented interface to XML
                                objects
    * XMLParameterListReader  - ParameterList input from XML
    * XMLParameterListWriter  - ParameterList output to XML
    * XMLInputSource          - Base class for converting a stream
                                to XML
    * FileInputSource         - Class for converting file contents
                                to XML
    * StringInputSource       - Class for converting string contents
                                to XML
    * ScalarTraits            - Function factory for ScalarTraits<...>
                                classes
    * Time                    - Wall-clock timer class

The ParameterList class matches string keys to arbitrarily-typed
values.  In python, the Teuchos.ParameterList is tightly integrated
with python dictionaries -- PyTrilinos methods that expect a
ParameterList will accept a python dictionary.
"
%enddef

%module(package      = "PyTrilinos",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %teuchos_docstring) Teuchos

// SWIG does not support wrapping nested classes.  We will %include
// the Teuchos::ParameterList class, which has nested classes.  To
// suppress the swig warning that would otherwise result, we use the
// following:
#pragma SWIG nowarn=312

// Includes
%{
// Configuration includes
#include "PyTrilinos_config.h"
#include "Teuchos_ConfigDefs.hpp"

// Teuchos includes
#include "Teuchos_FILEstream.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_NullIteratorTraits.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterListExceptions.hpp"
#include "Teuchos_Time.hpp"

// Local includes
#include "NumPyImporter.h"

// Namespace flattening
using std::string;
using Teuchos::RCP;
%}

// Namespace flattening
using std::string;
using Teuchos::RCP;

// Standard exception handling
%include "exception.i"

// Global swig features
%feature("autodoc", "1");
%feature("compactdefaultargs");

// Include Teuchos documentation
%include "Teuchos_dox.i"

// C++ STL support.  If the wrapped class uses standard template
// library containers, the following %include wraps the containers
// and makes certain conversions seamless, such as between std::string
// and python strings.
%include "stl.i"
%include "std_except.i"

// Teuchos utilizes very little of numpy.i, so some of the fragments
// do not get automatically instantiated.  This forces the issue.
%include "numpy.i"
%fragment("NumPy_Fragments");

// General exception handling
%exception
{
  try
  {
    $action
    if (PyErr_Occurred()) SWIG_fail;
  }
  catch(Teuchos::Exceptions::InvalidParameterType & e)
  {
    SWIG_exception(SWIG_TypeError, e.what());
  }
  catch(Teuchos::Exceptions::InvalidParameter & e)
  {
    PyErr_SetString(PyExc_KeyError, e.what());
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

// General ignore directives
%ignore *::operator=;
%ignore *::print;

//Teuchos imports
%import "Teuchos_TypeNameTraits.hpp"
%import "Teuchos_NullIteratorTraits.hpp"
%import "Teuchos_RCPDecl.hpp"

// Teuchos includes
%include "Teuchos_Traits.i"
%include "Teuchos_Comm.i"
%include "Teuchos_ParameterList.i"
%include "Teuchos_XML.i"

//////////////////////////////////////
// Teuchos::Teuchos_Version support //
//////////////////////////////////////
%include "Teuchos_Version.hpp"
%pythoncode
%{
__version__ = Teuchos_Version().split()[2]
%}

///////////////////////////
// Teuchos::Time support //
///////////////////////////
%include "Teuchos_Time.hpp"

////////////////////////////////////////////////////////////////////////

// Extend the %extend_smart_pointer macro to handle derived classes.
// This has been specialized for RCP, because the constructor for
// creating an RCP has an additional argument: a boolean false.
%define %extend_RCP(Type...)
%typemap(in, noblock=1) const SWIGTYPE & SMARTPOINTER
(void* argp = 0, int res = 0)
{
  res = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
  if (!SWIG_IsOK(res))
  {
    res = SWIG_ConvertPtr($input, &argp, $descriptor(Type*), %convertptr_flags);
    if (!SWIG_IsOK(res))
    {
      %argument_fail(res, "$type", $symname, $argnum);
    }
    if (!argp)
    {
      %argument_nullref("$type", $symname, $argnum);
    }
    $1 = new $*ltype ( %reinterpret_cast( argp, Type* ), false );
  }
  else
  {
    if (!argp) { %argument_nullref("$type", $symname, $argnum); }
    $1 = %reinterpret_cast(argp, $ltype);
  }
}

%typecheck(1200) const SWIGTYPE & SMARTPOINTER
{
  static void * argp = 0;
  $1 = SWIG_CheckState(SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags))
    ? 1 : 0;
  if (!$1)
    $1 = SWIG_CheckState(SWIG_ConvertPtr($input, &argp, $descriptor(Type *),
					 %convertptr_flags)) ? 1 : 0;
}

%typemap(freearg) const SWIGTYPE & SMARTPOINTER
{
  delete $1;
}

%typemap(out) SWIGTYPE SMARTPOINTER
{
  $result = SWIG_NewPointerObj((void*)$1.get(), $descriptor(Type*), 1);
}

%extend_smart_pointer(RCP< Type >)
%template()           RCP< Type >;

%enddef

// These typemap macros allow developers to generate typemaps for any
// classes that are wrapped in RCP<...> and used as function or method
// arguments.
%define %teuchos_rcp_typemaps(Type...)

%extend_RCP(      Type)
%extend_RCP(const Type)

%enddef

// Apply the RCP typemap macros to selected Teuchos classes
%teuchos_rcp_typemaps(Teuchos::Time)

// Provide specialized typemaps for RCP<Teuchos::ParameterList>
%typemap(in) RCP<Teuchos::ParameterList>
(Teuchos::ParameterList * params = NULL)
{
  int    res  = 0;
  void * argp = NULL;
  if (PyDict_Check($input))
  {
    params = Teuchos::pyDictToNewParameterList($input);
    if (!params)
    {
      PyErr_SetString(PyExc_TypeError,
		      "Python dictionary cannot be converted to ParameterList");
      SWIG_fail;
    }
  }
  else
  {
    res = SWIG_ConvertPtr($input, &argp, $descriptor(Teuchos::ParameterList*), 0);
    if (!SWIG_IsOK(res))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Argument $argnum cannot be converted to ParameterList");
      SWIG_fail;
    }
    params = reinterpret_cast< Teuchos::ParameterList * >(argp);
  }
  $1 = new Teuchos::RCP<Teuchos::ParameterList> (params, false);
}
%typemap(freearg) RCP<Teuchos::ParameterList>
{
  if (params$argnum) delete $1;
}
%apply RCP<Teuchos::ParameterList>
{
  RCP<Teuchos::ParameterList>&,
  const RCP<Teuchos::ParameterList>,
  const RCP<Teuchos::ParameterList>&,
  Teuchos::RCP<Teuchos::ParameterList>,
  Teuchos::RCP<Teuchos::ParameterList>&,
  const Teuchos::RCP<Teuchos::ParameterList>,
  const Teuchos::RCP<Teuchos::ParameterList>&
};

// Turn off the exception handling
%exception;
