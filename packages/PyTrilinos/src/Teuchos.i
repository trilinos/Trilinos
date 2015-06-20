// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
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
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %teuchos_docstring) Teuchos

// Includes
%{
// Configuration includes
#include "PyTrilinos_config.h"
#ifdef HAVE_INTTYPES_H
#undef HAVE_INTTYPES_H
#endif
#ifdef HAVE_STDINT_H
#undef HAVE_STDINT_H
#endif
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_DLLExportMacro.h"
#include "Teuchos_iostream_helpers.hpp"

// Teuchos includes
#include "Teuchos_Version.hpp"
#include "Teuchos_NullIteratorTraits.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterListExceptions.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_DataAccess.hpp"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"

// Namespace flattening
using std::string;
using Teuchos::RCP;
%}

// Configuration macros for SWIG
%include "Teuchos_config.h"
%include "Teuchos_ConfigDefs.hpp"
%include "Teuchos_DLLExportMacro.h"
%include "Teuchos_iostream_helpers.hpp"
%include "PyTrilinos_config.h"

// Namespace flattening
using std::string;

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
%feature("director:except")
{
  if ($error != NULL) {
    throw Swig::DirectorMethodException();
  }
}
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
  catch (Swig::DirectorException & e)
  {
    SWIG_fail;
  }
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

// General ignore directives
%ignore *::operator=;
%ignore *::print;

// Teuchos::RCP support.  If a class is ever passed to or from a
// function or method wrapped by a Teuchos::RCP<>, then it should be
// stored internally as a Teuchos::RCP<> as well.  This is
// accomplished by %include-ing Teuchos_RCP.i and calling the provided
// macro %teuchos_rcp() on the class.
%include "Teuchos_RCP.i"
%teuchos_rcp(std::basic_ostream)
%teuchos_rcp(std::ostream)
%teuchos_rcp(std::vector< int, std::allocator< int > >)
%teuchos_rcp(Teuchos::SerialDenseMatrix< int, double >)

// Teuchos imports
%import "Teuchos_TypeNameTraits.hpp"
%import "Teuchos_NullIteratorTraits.hpp"

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
%teuchos_rcp(Teuchos::Time)
%include "Teuchos_Time.hpp"

/////////////////////////////////
// Teuchos::DataAccess support //
/////////////////////////////////
%include "Teuchos_DataAccess.hpp"

// Turn off the exception handling
%exception;
