// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2014) Sandia Corporation
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

%define %tpetra_docstring
"
PyTrilinos.Tpetra is the python interface to the Trilinos linear
algebra services package Tpetra:

    http://trilinos.sandia.gov/packages/tpetra

The purpose of Tpetra is to provide fundamental linear algebra
services to the rest of Trilinos.  These services include parallel
decomposition and communication, vectors and multivectors, graphs,
operators, and dense and sparse matrices.
"
%enddef

%module(package   = "PyTrilinos",
	directors = "1",
	docstring = %tpetra_docstring) Tpetra

%{
// PyTrilinos includes
#include "PyTrilinos_config.h"
#include "PyTrilinos_PythonException.h"
#include "PyTrilinos_Teuchos_Util.h"
#include "PyTrilinos_NumPy_Util.hpp"

// Import the numpy interface
#define NO_IMPORT_ARRAY
#include "numpy_include.h"

// Teuchos includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_ArrayRCP.hpp"
using Teuchos::RCP;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::ArrayRCP;

// Tpetra includes
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_CombineMode.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
%}

// Global swig features
%feature("autodoc", "1");
%feature("compactdefaultargs");

// SWIG standard library includes
using std::string;
%include "stl.i"

// SWIG NumPy interface file
%include "numpy.i"
%pythoncode
{
import numpy
}

// PyTrilinos configuration support
%include "PyTrilinos_config.h"

// Teuchos support
%import "Teuchos.i"
%include "Teuchos_Array.i"

// Include the standard exception handlers
%include "exception.i"

// General exception handling
%exception
{
  try
  {
    $action
    if (PyErr_Occurred()) SWIG_fail;
  }
  catch(PyTrilinos::PythonException &e)
  {
    e.restore();
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

// General ignore directives
%ignore *::operator[];
%ignore *::operator++;
%ignore *::operator--;

// Include Tpetra documentation
%include "Tpetra_dox.i"

////////////////////////////
// Tpetra version support //
////////////////////////////
%include "Tpetra_Version.hpp"
%pythoncode
%{
__version__ = version()
%}

/////////////////////////////////////////////
// Tpetra enumerations and typedef support //
/////////////////////////////////////////////
// Use %import and forward declarations to prevent SWIG warnings when
// we %include "Tpetra_ConfigDefs.hpp"
%import "Teuchos_config.h"
%import "Teuchos_ConfigDefs.hpp"
%import "Teuchos_ENull.hpp"
%import "Teuchos_Array.hpp"
%import "Teuchos_ArrayViewDecl.hpp"
%import "Teuchos_ArrayRCPDecl.hpp"
%import "Teuchos_Tuple.hpp"
%import "Teuchos_PtrDecl.hpp"
%import "Teuchos_OrdinalTraits.hpp"
namespace Teuchos
{
template< class T > RCP< T > rcp(T* p, bool owns_mem = true);
template< class T > RCP< T > rcpFromRef(T& r);
template< class T2, class T1 > RCP< T2 > rcp_const_cast(const RCP< T1 >& p1);
}
%include "Tpetra_ConfigDefs.hpp"
%include "Tpetra_CombineMode.hpp"

////////////////////////
// Tpetra Map support //
////////////////////////
%extend Tpetra::Map
{
  Map(Tpetra::global_size_t numGlobalElements,
      GlobalOrdinal indexBase,
      const Teuchos::RCP< const Teuchos::Comm< int > > & comm,
      Tpetra::LocalGlobal lg=GloballyDistributed)
  {
    return new Tpetra::Map< LocalOrdinal, GlobalOrdinal >(numGlobalElements,
                                                          indexBase,
                                                          comm,
                                                          lg);
  }

  Map(Tpetra::global_size_t numGlobalElements,
      size_t numLocalElements,
      GlobalOrdinal indexBase,
      const Teuchos::RCP< const Teuchos::Comm< int > > & comm)
  {
    return new Tpetra::Map< LocalOrdinal, GlobalOrdinal >(numGlobalElements,
                                                          numLocalElements,
                                                          indexBase,
                                                          comm);
  }

  Map(Tpetra::global_size_t numGlobalElements,
      PyObject * elementList,
      GlobalOrdinal indexBase,
      const Teuchos::RCP< const Teuchos::Comm< int > > & comm)
  {
    int is_new = 0;
    int type_code = PyTrilinos::NumPy_TypeCode< GlobalOrdinal >();
    PyArrayObject * npArray =
      obj_to_array_contiguous_allow_conversion(elementList,
                                               type_code,
                                               &is_new);
    if (!npArray) throw PyTrilinos::PythonException();
    Teuchos::ArrayView< GlobalOrdinal > elementArray =
      Teuchos::arrayView( (GlobalOrdinal*) array_data(npArray),
                          array_size(npArray, 0));
    return new Tpetra::Map< LocalOrdinal, GlobalOrdinal >(numGlobalElements,
                                                          elementArray,
                                                          indexBase,
                                                          comm);
  }

  PyObject * getLocalElement(GlobalOrdinal globalIndex)
  {
    LocalOrdinal localIndex = self->getLocalElement(globalIndex);
    if (localIndex == Teuchos::OrdinalTraits< LocalOrdinal >::invalid())
      return Py_BuildValue("");
    return SWIG_From_long(static_cast< long >(localIndex));
  }

  PyObject * getGlobalElement(LocalOrdinal localIndex)
  {
    GlobalOrdinal globalIndex = self->getGlobalElement(localIndex);
    if (globalIndex == Teuchos::OrdinalTraits< GlobalOrdinal >::invalid())
      return Py_BuildValue("");
    return SWIG_From_long(static_cast< long >(globalIndex));
  }

  PyObject * getRemoteIndexList(PyObject * globalIndexes)
  {
    // Variable initialization
    int is_new = 0;
    int type_code = PyTrilinos::NumPy_TypeCode< GlobalOrdinal >();
    npy_intp dims[1];
    PyObject * globalArray = NULL;
    PyObject * nodeArray   = NULL;
    PyObject * localArray  = NULL;
    PyObject * resultObj   = PyTuple_New(3);
    Teuchos::ArrayView< const GlobalOrdinal > globalList;
    Teuchos::ArrayView< int >                 nodeList;
    Teuchos::ArrayView< LocalOrdinal >        localList;
    Tpetra::LookupStatus result;

    // Check the input argument
    if (!(is_array(globalIndexes) || PySequence_Check(globalIndexes)))
    {
      PyErr_SetString(PyExc_ValueError,
                      "Argument to Tpetra.Map.getRemoteIndexList() is "
                      "not a NumPy array or Python sequence");
      goto fail;
    }

    // Convert the input argument
    globalArray = (PyObject *)
      obj_to_array_contiguous_allow_conversion(globalIndexes,
                                               type_code,
                                               &is_new);
    if (!globalArray) goto fail;
    if (!require_dimensions((PyArrayObject *) globalArray,1)) goto fail;

    // Initialize the output NumPy arrays and tuple
    dims[0]    = array_size(globalArray, 0);
    type_code  = PyTrilinos::NumPy_TypeCode< LocalOrdinal >();
    nodeArray  = PyArray_SimpleNew(1, dims, NPY_INT);
    localArray = PyArray_SimpleNew(1, dims, type_code);
    PyTuple_SET_ITEM(resultObj, 0, nodeArray );
    PyTuple_SET_ITEM(resultObj, 1, localArray);

    // Initialize the Teuchos::ArrayViews
    globalList = Teuchos::ArrayView< const GlobalOrdinal >
      ((const GlobalOrdinal*) array_data(globalArray),
       (Teuchos::ArrayView< GlobalOrdinal >::size_type) dims[0]);
    nodeList = Teuchos::ArrayView< int >
      ((int*) array_data(nodeArray),
       (Teuchos::ArrayView< int >::size_type) dims[0]);
    localList = Teuchos::ArrayView< LocalOrdinal >
      ((LocalOrdinal*) array_data(localArray),
       (Teuchos::ArrayView< LocalOrdinal >::size_type) dims[0]);

    // Call the method
    result = self->getRemoteIndexList(globalList,
                                      nodeList,
                                      localList);

    // Final conversions for output object
    PyTuple_SET_ITEM(resultObj, 2, SWIG_From_long(static_cast< long >(result)));
    return resultObj;

    // Error handling
  fail:
    if (is_new) Py_XDECREF(globalArray);
    Py_XDECREF(nodeArray);
    Py_XDECREF(localArray);
    Py_XDECREF(resultObj);
    return NULL;
  }

  std::string __str__()
  {
    return self->description();
  }
}
%ignore Tpetra::Map::Map;
%ignore Tpetra::Map::getLocalElement;
%ignore Tpetra::Map::getGlobalElement;
%ignore Tpetra::Map::getRemoteIndexList;
%include "Tpetra_Map_decl.hpp"
%teuchos_rcp(Tpetra::Map< long, long >)
%template(Map_default) Tpetra::Map< long, long >;
%pythoncode
{
Map = Map_default
}

/////////////////////////////
// Tpetra Transfer support //
/////////////////////////////
%include "Tpetra_Details_Transfer.hpp"
%teuchos_rcp(Tpetra::Details::Transfer< long, long, KokkosClassic::DefaultNode::DefaultNodeType >)
%template(Transfer_default) Tpetra::Details::Transfer< long, long, KokkosClassic::DefaultNode::DefaultNodeType >;

///////////////////////////
// Tpetra Export support //
///////////////////////////
%include "Tpetra_Export_decl.hpp"
%teuchos_rcp(Tpetra::Export< long, long, KokkosClassic::DefaultNode::DefaultNodeType >)
%template(Export_default) Tpetra::Export< long, long, KokkosClassic::DefaultNode::DefaultNodeType >;
%pythoncode
{
Export = Export_default
}

///////////////////////////
// Tpetra Import support //
///////////////////////////
%include "Tpetra_Import_decl.hpp"
%teuchos_rcp(Tpetra::Import< long, long, KokkosClassic::DefaultNode::DefaultNodeType >)
%template(Import_default) Tpetra::Import< long, long, KokkosClassic::DefaultNode::DefaultNodeType >;
%pythoncode
{
Import = Import_default
}

//////////////////////////////////
// Tpetra SrcDistObject support //
//////////////////////////////////
%teuchos_rcp(Tpetra::SrcDistObject)
%include "Tpetra_SrcDistObject.hpp"

///////////////////////////////
// Tpetra DistObject support //
///////////////////////////////
%ignore Tpetra::removeEmptyProcessesInPlace;
%include "Tpetra_DistObject_decl.hpp"
%teuchos_rcp(Tpetra::DistObject< int, long, long >)
%template(DistObject_int) Tpetra::DistObject< int, long, long >;

////////////////////////////////
// Tpetra MultiVector support //
////////////////////////////////
// %include "Tpetra_MultiVector_decl.hpp"
// %teuchos_rcp(Tpetra::MultiVector< int, long, long >)
// %template(MultiVector_int) Tpetra::MultiVector< int, long, long >;

// Turn off exception handling
%exception;
