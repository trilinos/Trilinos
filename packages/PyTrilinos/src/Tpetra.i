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
#include "PyTrilinos_PythonException.hpp"
#include "PyTrilinos_Teuchos_Util.hpp"
#include "PyTrilinos_NumPy_Util.hpp"

// Import the numpy interface
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"

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

// Define a shortcut for the default Kokkos node
%inline
%{
  typedef KokkosClassic::DefaultNode::DefaultNodeType KokkosDefaultNode;
%}

////////////////////////////////////////////////////////////
// Tpetra configuration, enumerations and typedef support //
////////////////////////////////////////////////////////////
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
%include "KokkosCore_config.h"
%include "Kokkos_Macros.hpp"
%include "TpetraCore_config.h"
%include "Tpetra_ConfigDefs.hpp"
%include "Tpetra_CombineMode.hpp"

////////////////////////////
// Tpetra version support //
////////////////////////////
%include "Tpetra_Version.hpp"
%pythoncode
%{
__version__ = version()
%}

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
// If I do a %include "Tpetra_Details_Transfer.hpp", SWIG aborts with
// a trap 6.  This appears to be related to the template parameters of
// class Tpetra::Details::Transfer using absolute namespaces
// ::Tpetra::Map<>::... as default parameters.  I reproduce the
// Transfer class here without the absolute namespace defaults.  This
// could cause problems if the Transfer class changes.
// %include "Tpetra_Details_Transfer.hpp"
namespace Tpetra
{
class Distributor;
namespace Details
{
template <class LO = Tpetra::Map<>::local_ordinal_type,
          class GO = typename Tpetra::Map<LO>::global_ordinal_type,
          class NT = typename Tpetra::Map<LO, GO>::node_type>
class Transfer : public Teuchos::Describable
{
public:
  virtual ~Transfer () {}
  typedef ::Tpetra::Map<LO, GO, NT> map_type;
  virtual size_t getNumSameIDs() const = 0;
  virtual size_t getNumPermuteIDs () const = 0;
  virtual Teuchos::ArrayView<const LO> getPermuteFromLIDs () const = 0;
  virtual Teuchos::ArrayView<const LO> getPermuteToLIDs () const = 0;
  virtual size_t getNumRemoteIDs () const = 0;
  virtual Teuchos::ArrayView<const LO> getRemoteLIDs () const = 0;
  virtual size_t getNumExportIDs () const = 0;
  virtual Teuchos::ArrayView<const LO> getExportLIDs () const = 0;
  virtual Teuchos::ArrayView<const int> getExportPIDs () const = 0;
  virtual Teuchos::RCP<const map_type> getSourceMap () const = 0;
  virtual Teuchos::RCP<const map_type> getTargetMap () const = 0;
  virtual ::Tpetra::Distributor& getDistributor () const = 0;
};
} // namespace Details
} // namespace Tpetra
%teuchos_rcp(Tpetra::Details::Transfer< long, long, KokkosDefaultNode >)
%template(Transfer_default)
    Tpetra::Details::Transfer< long, long, KokkosDefaultNode >;

///////////////////////////
// Tpetra Export support //
///////////////////////////
%include "Tpetra_Export_decl.hpp"
%teuchos_rcp(Tpetra::Export< long, long, KokkosDefaultNode >)
%template(Export_default) Tpetra::Export< long, long, KokkosDefaultNode >;
%pythoncode
{
Export = Export_default
}

///////////////////////////
// Tpetra Import support //
///////////////////////////
%include "Tpetra_Import_decl.hpp"
%teuchos_rcp(Tpetra::Import< long, long, KokkosDefaultNode >)
%template(Import_default) Tpetra::Import< long, long, KokkosDefaultNode >;
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
%include "Tpetra_KokkosRefactor_DistObject_decl.hpp"

////////////////////////////////
// Tpetra MultiVector support //
////////////////////////////////
%ignore Tpetra::MultiVector::getLocalMV;
%typemap(in)
  (const Teuchos::ArrayView<const Scalar> & A, size_t LDA, size_t NumVectors)
  (PyArrayObject * array=NULL, int is_new=0)
{
  array =
    obj_to_array_contiguous_allow_conversion($input,
                                             PyTrilinos::NumPy_TypeCode<Scalar>(),
                                             &is_new);
  if (!array || !require_dimensions(array, 2))
    SWIG_fail;
  $2 = array_stride(array, 0);
  $3 = array_size(array, 1);
  $1 = Teuchos::arrayView((Scalar*)array_data(array), $2, $3)
}
%typemap(freearg)
  (const Teuchos::ArrayView<const Scalar> & A, size_t LDA, size_t NumVectors)
{
  if (is_new$argnum && array$argnum)
  {
    Py_DECREF(array$argnum);
  }
}
%feature("notabstract") Tpetra::MultiVector;
%include "Tpetra_MultiVector_decl.hpp"
%include "Tpetra_KokkosRefactor_MultiVector_decl.hpp"

///////////////////////////
// Tpetra Vector support //
///////////////////////////
// %ignore Tpetra::Vector::getLocalMV;
%feature("notabstract") Tpetra::Vector;
%include "Tpetra_Vector_decl.hpp"
%include "Tpetra_KokkosRefactor_Vector_decl.hpp"

/////////////////////////////////////
// Explicit template instantiation //
/////////////////////////////////////
%define %tpetra_class( CLASS, SCALAR, SCALAR_NAME )
    %warnfilter(315) Tpetra::CLASS< SCALAR, long, long, KokkosDefaultNode >;
    %teuchos_rcp(Tpetra::CLASS< SCALAR, long, long, KokkosDefaultNode >)
    %template(CLASS ## _ ## SCALAR_NAME)
        Tpetra::CLASS< SCALAR, long, long, KokkosDefaultNode, true >;
%enddef

%define %tpetra_scalars( SCALAR, SCALAR_NAME)
    %tpetra_class( DistObject , SCALAR, SCALAR_NAME )
    %tpetra_class( MultiVector, SCALAR, SCALAR_NAME )
    %tpetra_class( Vector     , SCALAR, SCALAR_NAME )
%enddef

//////////////////////////////////////////////
// Concrete scalar types for Tpetra classes //
//////////////////////////////////////////////
%tpetra_scalars(int   , int   )
%tpetra_scalars(long  , long  )
%tpetra_scalars(float , float )
%tpetra_scalars(double, double)

/////////////////////////////////////////////////////
// Python code that consolidates templated classes //
/////////////////////////////////////////////////////
%pythoncode
{
  def MultiVector(*args, **kwargs):
    dtype = kwargs.get("dtype", "int64")
    if type(dtype) == str:
      dtype = numpy.dtype(dtype)
    if dtype.type is numpy.int32:
      result = MultiVector_int(*args)
    elif dtype.type is numpy.int64:
      result = MultiVector_long(*args)
    elif dtype.type is numpy.float32:
      result = MultiVector_float(*args)
    elif dtype.type is numpy.flat64:
      result = MultiVector_double(*args)
    else:
      raise TypeError("Unsupported or unrecognized dtype = %s" %
                      str(dtype))
    return result

  def Vector(*args, **kwargs):
    dtype = kwargs.get("dtype", "int64")
    if type(dtype) == str:
      dtype = numpy.dtype(dtype)
    if dtype.type is numpy.int32:
      result = Vector_int(*args)
    elif dtype.type is numpy.int64:
      result = Vector_long(*args)
    elif dtype.type is numpy.float32:
      result = Vector_float(*args)
    elif dtype.type is numpy.flat64:
      result = Vector_double(*args)
    else:
      raise TypeError("Unsupported or unrecognized dtype = %s" %
                      str(dtype))
    return result
}

// Turn off exception handling
%exception;
