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
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_any.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_XMLObjectImplem.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_XMLInputSource.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_StringInputSource.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Time.hpp"

// Teuchos python interface includes
#include "Teuchos_PythonParameter.h"

// Namespace flattening
using std::string;
using Teuchos::RCP;
%}

// Namespace flattening
using std::string;
using Teuchos::RCP;

// Global swig features
%feature("autodoc", "1");
%feature("compactdefaultargs");

// Include Teuchos documentation
%include "Teuchos_dox.i"    // Doxygen-generated documentation
%include "Teuchos_doc.i"    // Manually written documentation

// C++ STL support.  If the wrapped class uses standard template
// library containers, the following %include wraps the containers
// and makes certain conversions seamless, such as between std::string
// and python strings.
%include "stl.i"
%include "std_except.i"

// Define macro for handling exceptions thrown by Teuchos methods and
// constructors
%define %teuchos_exception(className,methodName)
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
namespace Teuchos { class any; }
%import "Teuchos_TypeNameTraits.hpp"
%import "Teuchos_RCPDecl.hpp"
%import "Teuchos_ParameterEntry.hpp"
%import "Teuchos_XMLObjectImplem.hpp"
%import "Teuchos_PythonParameter.h"

//////////////////////////////////////
// Teuchos::Teuchos_Version support //
//////////////////////////////////////
%include "Teuchos_Version.hpp"
%pythoncode %{
__version__ = Teuchos_Version().split()[2]
%}

////////////////////////////////////
// Teuchos::ParameterList support //
////////////////////////////////////
%teuchos_exception(ParameterList, ParameterList)
%teuchos_exception(ParameterList, set)
%teuchos_exception(ParameterList, setParameters)
%teuchos_exception(ParameterList, get)
%teuchos_exception(ParameterList, sublist)
%teuchos_exception(ParameterList, type)
%teuchos_exception(ParameterList, __setitem__)
%teuchos_exception(ParameterList, update)
// There are a lot of extensions to the Teuchos::ParameterList class,
// so I put them all in their own file
%include "Teuchos_ParameterList_ext.i"
%ignore Teuchos::ParameterList::set;
%ignore Teuchos::ParameterList::setEntry;
%ignore Teuchos::ParameterList::get;
%ignore Teuchos::ParameterList::getPtr;
%ignore Teuchos::ParameterList::getEntryPtr;
%ignore Teuchos::ParameterList::sublist(const std::string &) const;
%ignore Teuchos::ParameterList::isType(const std::string &) const;
%ignore Teuchos::ParameterList::isType(const std::string &, any*) const;
%ignore Teuchos::ParameterList::unused(ostream &) const;
%ignore Teuchos::ParameterList::begin() const;
%ignore Teuchos::ParameterList::end() const;
%ignore Teuchos::ParameterList::entry(ConstIterator) const;
%ignore Teuchos::ParameterList::name(ConstIterator) const;
%include "Teuchos_ParameterList.hpp"
%typemap(in)
Teuchos::ParameterList &
(void *argp=0, int res=0, bool cleanup=false)
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
%typecheck(200)
Teuchos::ParameterList &
{
  // Accept PyDicts or ParameterLists
  void * argp = NULL;
  $1 = PyDict_Check($input) ? 1 : 0;
  if (!$1) if (SWIG_CheckState(SWIG_Python_ConvertPtr($input, &argp,
						      $1_descriptor, 0))) $1 = 1;
}
%typemap(freearg)
Teuchos::ParameterList &
{
  if (cleanup$argnum && $1) delete $1;
}

////////////////////////////////////////////
// Teuchos::ParameterListAcceptor support //
////////////////////////////////////////////
%ignore Teuchos::ParameterListAcceptor::setParameterList(RCP<Teuchos::ParameterList > const &);
%ignore Teuchos::ParameterListAcceptor::getParameterList();
%ignore Teuchos::ParameterListAcceptor::unsetParameterList();
%ignore Teuchos::ParameterListAcceptor::getParameterList() const;
%ignore Teuchos::ParameterListAcceptor::getValidParameters() const;
%extend Teuchos::ParameterListAcceptor {

  // The ParameterListAcceptor Class has the following virtual
  // functions with default implementations: getParameterList() const;
  // and getValidParameters() const.  These both return RCP<const
  // ParameterList> objects, which presents a problem: there are
  // typemaps for RCP< > and for ParameterList, but combinations of
  // the two must be handled in a "brute force" manner.

  const Teuchos::ParameterList * getParameterList() const {
    RCP<const Teuchos::ParameterList> p_plist = self->getParameterList();
    return p_plist.get();
  }

  const Teuchos::ParameterList * getValidParameters() const {
    RCP<const Teuchos::ParameterList> p_plist = self->getValidParameters();
    return p_plist.get();
  }
}
%include "Teuchos_ParameterListAcceptor.hpp"

////////////////////////////////
// Teuchos::XMLObject support //
////////////////////////////////
%teuchos_exception(XMLObject, deepCopy         )
%teuchos_exception(XMLObject, checkTag         )
%teuchos_exception(XMLObject, getTag           )
%teuchos_exception(XMLObject, getAttribute     )
%teuchos_exception(XMLObject, getContentLine   )
%teuchos_exception(XMLObject, getRequired      )
%teuchos_exception(XMLObject, getRequiredDouble)
%teuchos_exception(XMLObject, getRequiredInt   )
%teuchos_exception(XMLObject, getRequiredBool  )
%teuchos_exception(XMLObject, getChild         )
%teuchos_exception(XMLObject, addAttribute     )
%teuchos_exception(XMLObject, addContentLine   )
%teuchos_exception(XMLObject, addRequired      )
%teuchos_exception(XMLObject, addRequiredDouble)
%teuchos_exception(XMLObject, addRequiredInt   )
%teuchos_exception(XMLObject, addRequiredBool  )
%teuchos_exception(XMLObject, addChild         )
%teuchos_exception(XMLObject, toString         )
%ignore Teuchos::XMLObject::XMLObject(XMLObjectImplem*);
%extend Teuchos::XMLObject {
  std::string __str__() const {
    try {
      return self->toString();
    }
    catch(std::logic_error e) {
      return std::string("");
    }
  }
}
%include "Teuchos_XMLObject.hpp"

/////////////////////////////////////////////
// Teuchos::XMLParameterListReader support //
/////////////////////////////////////////////
%include "Teuchos_XMLParameterListReader.hpp"

/////////////////////////////////////////////
// Teuchos::XMLParameterListWriter support //
/////////////////////////////////////////////
%include "Teuchos_XMLParameterListWriter.hpp"

/////////////////////////////////////
// Teuchos::XMLInputSource support //
/////////////////////////////////////
%teuchos_exception(XMLInputSource, getObject)
%ignore Teuchos::XMLInputSource::stream() const;
%include "Teuchos_XMLInputSource.hpp"

//////////////////////////////////////
// Teuchos::FileInputSource support //
//////////////////////////////////////
%ignore Teuchos::FileInputSource::stream() const;
%include "Teuchos_FileInputSource.hpp"

////////////////////////////////////////
// Teuchos::StringInputSource support //
////////////////////////////////////////
%ignore Teuchos::StringInputSource::stream() const;
%include "Teuchos_StringInputSource.hpp"

///////////////////////////////////
// Teuchos::ScalarTraits support //
///////////////////////////////////
%include "Teuchos_ScalarTraits.hpp"
%template(ScalarTraitsFloat ) Teuchos::ScalarTraits< float  >;
%template(ScalarTraitsDouble) Teuchos::ScalarTraits< double >;
%pythoncode %{
def ScalarTraits(scalarType):
    """
    ScalarTraits(str scalarType) -> ScalarTraits<...>

    The scalarType argument is for specifying the type of scalar for
    which traits are requested.  Limited NumPy-style type
    specification is supported: 'f' for float and 'd' for double.
    """
    if scalarType == 'f': return ScalarTraitsFloat()
    if scalarType == 'd': return ScalarTraitsDouble()
    raise NotImplementedError, "ScalarTraits for " + repr(scalarType) + " not supported"
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

%typecheck(1200) const SWIGTYPE & SMARTPOINTER
{
  static void * argp = 0;
  $1 = SWIG_CheckState(SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags)) ? 1 : 0;
  if (!$1)
    $1 = SWIG_CheckState(SWIG_ConvertPtr($input, &argp, $descriptor(Type *),
					 %convertptr_flags)) ? 1 : 0;
}

%typemap(freearg) const SWIGTYPE & SMARTPOINTER
{
  delete $1;
}

//%typemap(out) SWIGTYPE SMARTPOINTER
//{
//  $result = SWIG_NewPointerObj($1.get(), $descriptor(Type*), 1);
//}

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
  int                      res          = 0;
  void                  *  argp         = NULL;
  static swig_type_info *  swig_TPL_ptr = SWIG_TypeQuery("Teuchos::ParameterList *");
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
      res = SWIG_ConvertPtr($input, &argp, swig_TPL_ptr, 0);
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
