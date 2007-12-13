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
#include "Teuchos_NullIteratorTraits.hpp"
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
#include "Teuchos_VerbosityLevel.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_LabeledObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ReductionOp.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#include "Teuchos_CommHelpers.hpp"

// Teuchos python interface includes
#include "Teuchos_PythonParameter.h"

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
%include "Teuchos_dox.i"    // Doxygen-generated documentation
%include "Teuchos_doc.i"    // Manually written documentation

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
  try {
    $action
    if (PyErr_Occurred()) SWIG_fail;
  } catch(Teuchos::Exceptions::InvalidParameterType & e) {
    SWIG_exception(SWIG_TypeError, e.what());
  } catch(Teuchos::Exceptions::InvalidParameter & e) {
    PyErr_SetString(PyExc_KeyError, e.what());
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch(...) {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

// General ignore directives
%ignore *::operator=;
%ignore *::print;

//Teuchos imports
namespace Teuchos { class any; }
%import "Teuchos_TypeNameTraits.hpp"
%import "Teuchos_NullIteratorTraits.hpp"
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

/////////////////////////////////////
// Teuchos::VerbosityLevel support //
/////////////////////////////////////
%rename(verbosityLevelToString) Teuchos::toString;
%include "Teuchos_VerbosityLevel.hpp"

///////////////////////////////////
// Teuchos::FancyOStream support //
///////////////////////////////////
%include "Teuchos_FancyOStream.hpp"

////////////////////////////////////
// Teuchos::LabeledObject support //
////////////////////////////////////
%include "Teuchos_LabeledObject.hpp"

//////////////////////////////////
// Teuchos::Describable support //
//////////////////////////////////
%include "Teuchos_Describable.hpp"

//////////////////////////////////
// Teuchos::ReductionOp support //
//////////////////////////////////
%include "Teuchos_ReductionOp.hpp"

///////////////////////////
// Teuchos::Comm support //
///////////////////////////
%extend Teuchos::Comm
{
  PyObject * broadcast(int rootRank, PyObject * buffer) const
  {
    PyArrayObject * bufferArray = obj_to_array_no_conversion(buffer, NPY_NOTYPE);
    if (!bufferArray || !require_contiguous(bufferArray)) return NULL;
    Ordinal bytes = static_cast<Ordinal>(PyArray_NBYTES(bufferArray));
    char * bufferVals = (char*) array_data(bufferArray);
    self->broadcast(rootRank, bytes, bufferVals);
    return Py_BuildValue("");
  }
  PyObject * gatherAll(PyObject * sendBuffer) const
  {
    int     is_new_object  = 0;
    Ordinal sendBytes      = 0;
    Ordinal recvBytes      = 0;
    int     recvNd         = 0;
    int     type           = 0;
    char *  sendBufferVals = NULL;
    char *  recvBufferVals = NULL;
    PyArrayObject * recvBufferArray = NULL;
    PyArrayObject * sendBufferArray =
      obj_to_array_contiguous_allow_conversion(sendBuffer, NPY_NOTYPE, &is_new_object);
    if (!sendBufferArray) goto fail;
    sendBytes = static_cast<Ordinal>(PyArray_NBYTES(sendBufferArray));
    recvBytes = sendBytes * self->getSize();
    recvNd    = array_numdims(sendBufferArray) + 1;
    type      = array_type(sendBufferArray);
    { // Scope this to make recvDims temporary
      npy_intp * recvDims = new npy_intp[recvNd];
      recvDims[0] = self->getSize();
      for (int i=1; i<recvNd; ++i) recvDims[i] = array_size(sendBufferArray,i-1);
      recvBufferArray = (PyArrayObject*) PyArray_SimpleNew(recvNd, recvDims, type);
      if (!recvBufferArray) goto fail;
      delete [] recvDims;
    }
    sendBufferVals = (char*) array_data(sendBufferArray);
    recvBufferVals = (char*) array_data(recvBufferArray);
    self->gatherAll(sendBytes, sendBufferVals, recvBytes, recvBufferVals);
    if (is_new_object && sendBufferArray) Py_DECREF(sendBufferArray);
    return PyArray_Return(recvBufferArray);
  fail:
    if (is_new_object && sendBufferArray) Py_DECREF(sendBufferArray);
    return NULL;
  }
  PyObject * reduceAll(Teuchos::EReductionType reductOp, PyObject * sendBuffer) const
  {
    int     is_new_object = 0;
    Ordinal bytes         = 0;
    int     type          = 0;
    char * sendBufferVals    = NULL;
    char * globalReductsVals = NULL;
    PyArrayObject * globalReductsArray = NULL;
    PyArrayObject * sendBufferArray    =
      obj_to_array_contiguous_allow_conversion(sendBuffer, NPY_NOTYPE, &is_new_object);
    if (!sendBufferArray) goto fail;
    bytes = static_cast<Ordinal>(PyArray_NBYTES(sendBufferArray));
    type  = array_type(sendBufferArray);
    globalReductsArray = (PyArrayObject*)
      PyArray_SimpleNew(array_numdims(sendBufferArray),
			array_dimensions(sendBufferArray), type);
    PyArray_FILLWBYTE(globalReductsArray, 0);
    sendBufferVals    = (char*) array_data(sendBufferArray);
    globalReductsVals = (char*) array_data(globalReductsArray);
    Teuchos::reduceAll(*self, reductOp, bytes, sendBufferVals, globalReductsVals);
    if (is_new_object && sendBufferArray) Py_DECREF(sendBufferArray);
    return PyArray_Return(globalReductsArray);
  fail:
    if (is_new_object && sendBufferArray) Py_DECREF(sendBufferArray);
    return NULL;
  }
  PyObject * scan(Teuchos::EReductionType reductOp, PyObject * sendBuffer) const
  {
    int     is_new_object = 0;
    Ordinal bytes         = 0;
    int     type          = 0;
    char * sendBufferVals  = NULL;
    char * scanReductsVals = NULL;
    PyArrayObject * scanReductsArray = NULL;
    PyArrayObject * sendBufferArray  =
      obj_to_array_contiguous_allow_conversion(sendBuffer, NPY_NOTYPE, &is_new_object);
    if (!sendBufferArray) goto fail;
    bytes = static_cast<Ordinal>(PyArray_NBYTES(sendBufferArray));
    type  = array_type(sendBufferArray);
    scanReductsArray = (PyArrayObject*)
      PyArray_SimpleNew(array_numdims(sendBufferArray),
			array_dimensions(sendBufferArray), type);
    PyArray_FILLWBYTE(scanReductsArray, 0);
    sendBufferVals    = (char*) array_data(sendBufferArray);
    scanReductsVals = (char*) array_data(scanReductsArray);
    Teuchos::scan(*self, reductOp, bytes, sendBufferVals, scanReductsVals);
    if (is_new_object && sendBufferArray) Py_DECREF(sendBufferArray);
    return PyArray_Return(scanReductsArray);
  fail:
    if (is_new_object && sendBufferArray) Py_DECREF(sendBufferArray);
    return NULL;
  }
  std::string __str__() const
  {
    return self->description();
  }
}
%ignore Teuchos::Comm::broadcast;
%ignore Teuchos::Comm::gatherAll;
%ignore Teuchos::Comm::reduceAll;
%ignore Teuchos::Comm::scan;
%include "Teuchos_Comm.hpp"
%template(Comm_long) Teuchos::Comm<long>;
%pythoncode %{
Comm = Comm_long
%}

/////////////////////////////////
// Teuchos::SerialComm support //
/////////////////////////////////
%ignore Teuchos::SerialComm::broadcast;
%ignore Teuchos::SerialComm::gatherAll;
%ignore Teuchos::SerialComm::reduceAll;
%ignore Teuchos::SerialComm::scan;
%include "Teuchos_DefaultSerialComm.hpp"
%template(SerialComm_long) Teuchos::SerialComm<long>;
%pythoncode %{
SerialComm = SerialComm_long
PyComm = SerialComm
%}

//////////////////////////////////
// Teuchos::Comm Helper support //
//////////////////////////////////
%rename(reductionTypeToString) Teuchos::toString;
%include "Teuchos_CommHelpers.hpp"
%template(rank_long   ) Teuchos::rank<long>;
%template(size_long   ) Teuchos::size<long>;
%template(barrier_long) Teuchos::barrier<long>;
%pythoncode %{
rank    = rank_long
size    = size_long
barrier = barrier_long

def broadcast(comm, rootRank, buffer):
  """
  broadcast(Teuchos.Comm comm, int rootRank, numpy.ndarray buffer)

  Broadcast the contents of buffer from processor rootRank to all of the other
  processors.
  """
  comm.broadcast(rootRank, buffer)

def gatherAll(comm, buffer):
  """
  gatherAll(Teuchos.Comm comm, buffer) -> numpy.ndarray

  Gather the contents of buffer to all of the processors.
  """
  return comm.gatherAll(buffer)

def reduceAll(comm, reductOp, buffer):
  """
  Reduce the contents of buffer according to the operation designated by
  reductOp on all of the processors.
  """
  return comm.reduceAll(reductOp, buffer)

def scan(comm, reductOp, buffer):
  """
  Return the scan of the contents of buffer according to the operation
  designated by reductOp on each of the processors.
  """
  return comm.scan(reductOp, buffer)
%}

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

// Turn off the exception handling
%exception;
