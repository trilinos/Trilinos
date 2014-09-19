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

// Import the numpy interface
#define NO_IMPORT_ARRAY
#include "numpy_include.h"

// Teuchos includes
#include "Teuchos_RCP.hpp"
using Teuchos::RCP;

// Tpetra includes
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_Map.hpp"
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
//%include "Tpetra_dox.i"

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
    PyArrayObject * npArray =
      obj_to_array_contiguous_allow_conversion(elementList, NPY_LONG, &is_new);
    if (!npArray)
    {
      if (PyErr_Occurred()) throw PyTrilinos::PythonException();
      std::cout << "*** npArray is NULL, but no Python Exception raised ! ***" << std::endl;
    }
    Teuchos::ArrayView< GlobalOrdinal > elementArray =
      Teuchos::arrayView( (GlobalOrdinal*) array_data(npArray),
                          array_size(npArray, 0));
    return new Tpetra::Map< LocalOrdinal, GlobalOrdinal >(numGlobalElements,
                                                          elementArray,
                                                          indexBase,
                                                          comm);
  }
}
%ignore Tpetra::Map::Map;
%include "Tpetra_Map_decl.hpp"
%teuchos_rcp(Tpetra::Map< long, long >)
%template(Map_default) Tpetra::Map< long, long >;
%pythoncode
{
Map = Map_default
}

// Turn off exception handling
%exception;
