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
// PyTrilinos include files
#include "PyTrilinos_config.h"
#include "PyTrilinos_PythonException.hpp"
#include "PyTrilinos_NumPy_Util.hpp"
#include "PyTrilinos_Teuchos_Util.hpp"
#include "PyTrilinos_DAP.hpp"

// Import the numpy interface
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"

// Teuchos include files
#include "PyTrilinos_Teuchos_Headers.hpp"
using Teuchos::RCP;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::ArrayRCP;

// Tpetra include files
#include "PyTrilinos_Tpetra_Headers.hpp"

#ifdef HAVE_DOMI
// Domi include files
#include "PyTrilinos_Domi_Headers.hpp"
#endif

%}

// Define shortcuts for the default Tpetra template types
%inline
%{
  typedef Tpetra::Details::DefaultTypes::scalar_type         DefaultScalarType;
  typedef Tpetra::Details::DefaultTypes::local_ordinal_type  DefaultLOType;
  typedef Tpetra::Details::DefaultTypes::global_ordinal_type DefaultGOType;
  typedef Tpetra::Details::DefaultTypes::node_type           DefaultNodeType;

%}

%{
namespace PyTrilinos
{

// Attempt to convert a PyObject to an RCP of a Tpetra::MultiVector.
// The input PyObject could be a wrapped Tpetra::MultiVector, or a
// wrapped Domi::MDVector, or an object that supports the DistArray
// Protocol, or, if the environment is serial, a simple NumPy array.
template< class Scalar >
Teuchos::RCP< Tpetra::MultiVector< Scalar,
                                   PYTRILINOS_LOCAL_ORD,
                                   PYTRILINOS_GLOBAL_ORD,
                                   DefaultNodeType > > *
convertPythonToTpetraMultiVector(PyObject * pyobj,
                                 int * newmem)
{
  // SWIG initialization
  static swig_type_info * swig_TMV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Tpetra::MultiVector< Scalar,PYTRILINOS_LOCAL_ORD,PYTRILINOS_GLOBAL_ORD,DefaultNodeType > >*");
  static swig_type_info * swig_DMDV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Domi::MDVector< Scalar,Domi::DefaultNode::DefaultNodeType > >*");
  //
  // Get the default communicator
  const Teuchos::RCP< const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();
  //
  // Result objects
  void *argp = 0;
  PyObject * distarray = 0;
  Teuchos::RCP< Tpetra::MultiVector< Scalar,
                                     PYTRILINOS_LOCAL_ORD,
                                     PYTRILINOS_GLOBAL_ORD,
                                     DefaultNodeType > > smartresult;
  Teuchos::RCP< Tpetra::MultiVector< Scalar,
                                     PYTRILINOS_LOCAL_ORD,
                                     PYTRILINOS_GLOBAL_ORD,
                                     DefaultNodeType > > * result;
#ifdef HAVE_DOMI
  Teuchos::RCP< Domi::MDVector< Scalar > > dmdv_rcp;
#endif
  *newmem = 0;
  //
  // Check if the Python object is a wrapped Tpetra::MultiVector
  int res = SWIG_ConvertPtrAndOwn(pyobj, &argp, swig_TMV_ptr, 0, newmem);
  if (SWIG_IsOK(res))
  {
    result =
      reinterpret_cast< Teuchos::RCP< Tpetra::MultiVector< Scalar,
                                                           PYTRILINOS_LOCAL_ORD,
                                                           PYTRILINOS_GLOBAL_ORD,
                                                           DefaultNodeType > > * >(argp);
    return result;
  }

#ifdef HAVE_DOMI
  //
  // Check if the Python object is a wrapped Domi::MDVector< Scalar >
  *newmem = 0;
  res = SWIG_ConvertPtrAndOwn(pyobj, &argp, swig_DMDV_ptr, 0, newmem);
  if (SWIG_IsOK(res))
  {
    dmdv_rcp =
      *reinterpret_cast< Teuchos::RCP< Domi::MDVector< Scalar > > * >(argp);
    try
    {
      smartresult =
        dmdv_rcp->template getTpetraMultiVectorView< PYTRILINOS_LOCAL_ORD,
                                                     PYTRILINOS_GLOBAL_ORD >();
      *newmem = *newmem | SWIG_CAST_NEW_MEMORY;
    }
    catch (Domi::TypeError & e)
    {
      PyErr_SetString(PyExc_TypeError, e.what());
      return NULL;
    }
    catch (Domi::MDMapNoncontiguousError & e)
    {
      PyErr_SetString(PyExc_ValueError, e.what());
      return NULL;
    }
    catch (Domi::MapOrdinalError & e)
    {
      PyErr_SetString(PyExc_IndexError, e.what());
      return NULL;
    }
    result = new Teuchos::RCP< Tpetra::MultiVector< Scalar,
                                                    PYTRILINOS_LOCAL_ORD,
                                                    PYTRILINOS_GLOBAL_ORD,
                                                    DefaultNodeType > >(smartresult);
    return result;
  }
  //
  // Check if the Python object supports the DistArray Protocol
  if (PyObject_HasAttrString(pyobj, "__distarray__"))
  {
    try
    {
      if (!(distarray = PyObject_CallMethod(pyobj, (char*) "__distarray__", (char*) "")))
        return NULL;
      DistArrayProtocol dap(distarray);
      dmdv_rcp = convertToMDVector< Scalar >(comm, dap);
      Py_DECREF(distarray);
    }
    catch (PythonException & e)
    {
      e.restore();
      return NULL;
    }
    try
    {
      smartresult =
        dmdv_rcp->template getTpetraMultiVectorView< PYTRILINOS_LOCAL_ORD,
                                                     PYTRILINOS_GLOBAL_ORD>();
      *newmem = SWIG_CAST_NEW_MEMORY;
    }
    catch (Domi::TypeError & e)
    {
      PyErr_SetString(PyExc_TypeError, e.what());
      return NULL;
    }
    catch (Domi::MDMapNoncontiguousError & e)
    {
      PyErr_SetString(PyExc_ValueError, e.what());
      return NULL;
    }
    catch (Domi::MapOrdinalError & e)
    {
      PyErr_SetString(PyExc_IndexError, e.what());
      return NULL;
    }
    result =
      new Teuchos::RCP< Tpetra::MultiVector< Scalar,
                                             PYTRILINOS_LOCAL_ORD,
                                             PYTRILINOS_GLOBAL_ORD,
                                             DefaultNodeType > >(smartresult);
    return result;
  }
#endif

  //
  // Check if the environment is serial, and if so, check if the
  // Python object is a NumPy array
  if (comm->getSize() == 1)
  {
    if (PyArray_Check(pyobj))
    {
      PyArrayObject * array =
        (PyArrayObject*) PyArray_ContiguousFromObject(pyobj, NumPy_TypeCode< Scalar >(), 0, 0);
      if (!array) return NULL;
      size_t numVec, vecLen;
      int ndim = PyArray_NDIM(array);
      if (ndim == 1)
      {
        numVec = 1;
        vecLen = PyArray_DIM(array, 0);
      }
      else
      {
        numVec = PyArray_DIM(array, 0);
        vecLen = 1;
        for (int i=1; i < ndim; ++i) vecLen *= PyArray_DIM(array, i);
      }
      Scalar * data = (Scalar*) PyArray_DATA(array);
      Teuchos::ArrayView< Scalar > arrayView(data, vecLen*numVec);
      Teuchos::RCP< const Tpetra::Map< PYTRILINOS_LOCAL_ORD,
                                       PYTRILINOS_GLOBAL_ORD,DefaultNodeType > >
        map = Teuchos::rcp(new Tpetra::Map< PYTRILINOS_LOCAL_ORD,
                                            PYTRILINOS_GLOBAL_ORD,
                                            DefaultNodeType >(vecLen, 0, comm));
      smartresult =
        Teuchos::rcp(new Tpetra::MultiVector< Scalar,
                                              PYTRILINOS_LOCAL_ORD,
                                              PYTRILINOS_GLOBAL_ORD,
                                              DefaultNodeType >(map,
                                                                arrayView,
                                                                vecLen,
                                                                numVec));
      result =
        new Teuchos::RCP< Tpetra::MultiVector< Scalar,
                                               PYTRILINOS_LOCAL_ORD,
                                               PYTRILINOS_GLOBAL_ORD,
                                               DefaultNodeType > >(smartresult);
      return result;
    }
  }
  //
  // If we get to this point, then none of our known converters will
  // work, so it is time to set a Python error
  PyErr_Format(PyExc_TypeError, "Could not convert argument of type '%s'\n"
               "to a Tpetra::MultiVector",
               convertPyStringToChar(PyObject_Str(PyObject_Type(pyobj))));
  return NULL;
}

////////////////////////////////////////////////////////////////////////

// Attempt to convert a PyObject to an RCP to a Tpetra::Vector.  The
// input PyObject could be a wrapped Tpetra::Vector, or a wrapped
// Domi::MDVector, or an object that supports the DistArray Protocol,
// or, if the environment is serial, a simple NumPy array.
template< class Scalar >
Teuchos::RCP< Tpetra::Vector< Scalar,
                              PYTRILINOS_LOCAL_ORD,
                              PYTRILINOS_GLOBAL_ORD,
                              DefaultNodeType > > *
convertPythonToTpetraVector(PyObject * pyobj,
                            int * newmem)
{
  // SWIG initialization
  static swig_type_info * swig_TV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Tpetra::Vector< Scalar,PYTRILINOS_LOCAL_ORD,PYTRILINOS_GLOBAL_ORD,DefaultNodeType > >*");
  static swig_type_info * swig_DMDV_ptr =
    SWIG_TypeQuery("Teuchos::RCP< Domi::MDVector< Scalar,Domi::DefaultNode::DefaultNodeType > >*");
  //
  // Get the default communicator
  const Teuchos::RCP< const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();
  //
  // Result objects
  void *argp = 0;
  PyObject * distarray = 0;
  Teuchos::RCP< Tpetra::Vector< Scalar,
                                PYTRILINOS_LOCAL_ORD,
                                PYTRILINOS_GLOBAL_ORD,
                                DefaultNodeType > > smartresult;
  Teuchos::RCP< Tpetra::Vector< Scalar,
                                PYTRILINOS_LOCAL_ORD,
                                PYTRILINOS_GLOBAL_ORD,
                                DefaultNodeType > > * result;
#ifdef HAVE_DOMI
  Teuchos::RCP< Domi::MDVector< Scalar > > dmdv_rcp;
#endif
  *newmem = 0;
  //
  // Check if the Python object is a wrapped Tpetra::Vector
  int res = SWIG_ConvertPtrAndOwn(pyobj, &argp, swig_TV_ptr, 0, newmem);
  if (SWIG_IsOK(res))
  {
    result =
      reinterpret_cast< Teuchos::RCP< Tpetra::Vector< Scalar,
                                                      PYTRILINOS_LOCAL_ORD,
                                                      PYTRILINOS_GLOBAL_ORD,
                                                      DefaultNodeType > > * >(argp);
    return result;
  }

#ifdef HAVE_DOMI
  //
  // Check if the Python object is a wrapped Domi::MDVector< Scalar >
  *newmem = 0;
  res = SWIG_ConvertPtrAndOwn(pyobj, &argp, swig_DMDV_ptr, 0, newmem);
  if (SWIG_IsOK(res))
  {
    dmdv_rcp =
      *reinterpret_cast< Teuchos::RCP< Domi::MDVector< Scalar > > * >(argp);
    try
    {
      smartresult =
        dmdv_rcp->template getTpetraVectorView< PYTRILINOS_LOCAL_ORD,
                                                PYTRILINOS_GLOBAL_ORD >();
      *newmem = *newmem | SWIG_CAST_NEW_MEMORY;
    }
    catch (Domi::TypeError & e)
    {
      PyErr_SetString(PyExc_TypeError, e.what());
      return NULL;
    }
    catch (Domi::MDMapNoncontiguousError & e)
    {
      PyErr_SetString(PyExc_ValueError, e.what());
      return NULL;
    }
    catch (Domi::MapOrdinalError & e)
    {
      PyErr_SetString(PyExc_IndexError, e.what());
      return NULL;
    }
    result = new Teuchos::RCP< Tpetra::Vector< Scalar,
                                               PYTRILINOS_LOCAL_ORD,
                                               PYTRILINOS_GLOBAL_ORD,
                                               DefaultNodeType > >(smartresult);
    return result;
  }
  //
  // Check if the Python object supports the DistArray Protocol
  if (PyObject_HasAttrString(pyobj, "__distarray__"))
  {
    try
    {
      if (!(distarray = PyObject_CallMethod(pyobj, (char*) "__distarray__", (char*) "")))
        return NULL;
      DistArrayProtocol dap(distarray);
      dmdv_rcp = convertToMDVector< Scalar >(comm, dap);
      Py_DECREF(distarray);
    }
    catch (PythonException & e)
    {
      e.restore();
      return NULL;
    }
    try
    {
      smartresult =
        dmdv_rcp->template getTpetraVectorView< PYTRILINOS_LOCAL_ORD,
                                                PYTRILINOS_GLOBAL_ORD >();
      *newmem = SWIG_CAST_NEW_MEMORY;
    }
    catch (Domi::TypeError & e)
    {
      PyErr_SetString(PyExc_TypeError, e.what());
      return NULL;
    }
    catch (Domi::MDMapNoncontiguousError & e)
    {
      PyErr_SetString(PyExc_ValueError, e.what());
      return NULL;
    }
    catch (Domi::MapOrdinalError & e)
    {
      PyErr_SetString(PyExc_IndexError, e.what());
      return NULL;
    }
    result = new Teuchos::RCP< Tpetra::Vector< Scalar,
                                               PYTRILINOS_LOCAL_ORD,
                                               PYTRILINOS_GLOBAL_ORD,
                                               DefaultNodeType > >(smartresult);
    return result;
  }
#endif

  //
  // Check if the environment is serial, and if so, check if the
  // Python object is a NumPy array
  if (comm->getSize() == 1)
  {
    if (PyArray_Check(pyobj))
    {
      PyArrayObject * array =
        (PyArrayObject*) PyArray_ContiguousFromObject(pyobj, NumPy_TypeCode< Scalar >(), 0, 0);
      if (!array) return NULL;
      int ndim = PyArray_NDIM(array);
      size_t vecLen = 1;
      for (int i=1; i < ndim; ++i) vecLen *= PyArray_DIM(array, i);
      Scalar * data = (Scalar*) PyArray_DATA(array);
      Teuchos::ArrayView< const Scalar > arrayView(data, vecLen);
      Teuchos::RCP< const Tpetra::Map< PYTRILINOS_LOCAL_ORD,
                                       PYTRILINOS_GLOBAL_ORD,
                                       DefaultNodeType > > map =
        Teuchos::rcp(new Tpetra::Map< PYTRILINOS_LOCAL_ORD,
                                      PYTRILINOS_GLOBAL_ORD,
                                      DefaultNodeType >(vecLen, 0, comm));
      smartresult =
        Teuchos::rcp(new Tpetra::Vector< Scalar,
                                         PYTRILINOS_LOCAL_ORD,
                                         PYTRILINOS_GLOBAL_ORD,
                                         DefaultNodeType >(map, arrayView));
      result =
        new Teuchos::RCP< Tpetra::Vector< Scalar,
                                          PYTRILINOS_LOCAL_ORD,
                                          PYTRILINOS_GLOBAL_ORD,
                                          DefaultNodeType > >(smartresult);
      return result;
    }
  }
  //
  // If we get to this point, then none of our known converters will
  // work, so it is time to set a Python error
  PyErr_Format(PyExc_TypeError, "Could not convert argument of type '%s'\n"
               "to a Tpetra::Vector",
               convertPyStringToChar(PyObject_Str(PyObject_Type(pyobj))));
  return NULL;
}

}    // Namespace PyTrilinos

%}

// Global swig features
%feature("autodoc", "1");

// SWIG standard library include files
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

// Include Tpetra documentation
%include "Tpetra_dox.i"

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

//////////////////////////////
// Python utility functions //
//////////////////////////////
%pythoncode
%{
  def class_array_inplace_op(self, op_str, other):
    in_op = getattr(self.array, "__i"+op_str+"__")
    in_op(other.array)
    return self

  def class_array_math_op(self, op_str, other):
    # Initialize the result by calling the copy constructor
    result = self.__class__(self)
    # Get the equivalent in-place operator for the result
    in_op = getattr(result.array, "__i"+op_str+"__")
    try:
      in_op(other.array)
    except AttributeError:
      in_op(other)
    return result

  def class_array_rmath_op(self, op_str, other):
    # Initialize the result by calling the copy constructor
    result = self.__class__(self)
    indices = (slice(None),) * len(self.array.shape)
    result.array[indices] = other
    in_op = getattr(result.array, "__i"+op_str+"__")
    in_op(self.array)
    return result

  def class_array_add_math_ops(cls, op_str):
    setattr(cls,
            "__i"+op_str+"__",
            lambda self, other: class_array_inplace_op(self, op_str, other))
    setattr(cls,
            "__"+op_str+"__",
            lambda self, other: class_array_math_op(self, op_str, other))
    setattr(cls,
            "__r"+op_str+"__",
            lambda self, other: class_array_rmath_op(self, op_str, other))

  def class_array_add_math(cls):
    class_array_add_math_ops(cls, "add")
    class_array_add_math_ops(cls, "sub")
    class_array_add_math_ops(cls, "mul")
    class_array_add_math_ops(cls, "add")

  def class_array_comp_op(self, op_str, other):
    comp_op = getattr(self.array, "__"+op_str+"__")
    try:
      return comp_op(other.array)
    except AttributeError:
      return comp_op(other)

  def class_array_add_comp_op(cls, op_str):
    setattr(cls,
            "__"+op_str+"__",
            lambda self, other: class_array_comp_op(self, op_str, other))

  def class_array_add_comp(cls):
    class_array_add_comp_op(cls, "lt")
    class_array_add_comp_op(cls, "le")
    class_array_add_comp_op(cls, "eq")
    class_array_add_comp_op(cls, "ne")
    class_array_add_comp_op(cls, "gt")
    class_array_add_comp_op(cls, "ge")

%}

// General ignore directives
%ignore *::operator[];
%ignore *::operator++;
%ignore *::operator--;

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
%include "Kokkos_Macros.hpp"
%ignore KokkosClassic::ESweepDirection;
%include "Kokkos_ConfigDefs.hpp"
%include "Kokkos_DefaultNode.cpp"
%include "TpetraCore_config.h"
%include "TpetraClassic_config.h"
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

%extend Tpetra::Map< PYTRILINOS_LOCAL_ORD, PYTRILINOS_GLOBAL_ORD, DefaultNodeType >
{
  Map()
  {
    return new Tpetra::Map< PYTRILINOS_LOCAL_ORD,
                            PYTRILINOS_GLOBAL_ORD,
                            DefaultNodeType       >();
  }

  Map(Tpetra::global_size_t numGlobalElements,
      PYTRILINOS_GLOBAL_ORD indexBase,
      const Teuchos::RCP< const Teuchos::Comm< int > > & comm,
      Tpetra::LocalGlobal lg=GloballyDistributed)
  {
    return new Tpetra::Map< PYTRILINOS_LOCAL_ORD,
                            PYTRILINOS_GLOBAL_ORD,
                            DefaultNodeType       >(numGlobalElements,
                                                    indexBase,
                                                    comm,
                                                    lg);
  }

  Map(Tpetra::global_size_t numGlobalElements,
      size_t numLocalElements,
      PYTRILINOS_GLOBAL_ORD indexBase,
      const Teuchos::RCP< const Teuchos::Comm< int > > & comm)
  {
    return new Tpetra::Map< PYTRILINOS_LOCAL_ORD,
                            PYTRILINOS_GLOBAL_ORD,
                            DefaultNodeType       >(numGlobalElements,
                                                    numLocalElements,
                                                    indexBase,
                                                    comm);
  }

  Map(Tpetra::global_size_t numGlobalElements,
      PyObject * elementList,
      PYTRILINOS_GLOBAL_ORD indexBase,
      const Teuchos::RCP< const Teuchos::Comm< int > > & comm)
  {
    int is_new = 0;
    int type_code = PyTrilinos::NumPy_TypeCode< PYTRILINOS_GLOBAL_ORD >();
    PyArrayObject * npArray =
      obj_to_array_contiguous_allow_conversion(elementList,
                                               type_code,
                                               &is_new);
    if (!npArray) throw PyTrilinos::PythonException();
    Teuchos::ArrayView< PYTRILINOS_GLOBAL_ORD > elementArray =
      Teuchos::arrayView( (PYTRILINOS_GLOBAL_ORD*) array_data(npArray),
                          array_size(npArray, 0));
    return new Tpetra::Map< PYTRILINOS_LOCAL_ORD,
                            PYTRILINOS_GLOBAL_ORD,
                            DefaultNodeType       >(numGlobalElements,
                                                    elementArray,
                                                    indexBase,
                                                    comm);
  }

  PyObject * getLocalElement(PYTRILINOS_GLOBAL_ORD globalIndex)
  {
    PYTRILINOS_LOCAL_ORD localIndex = self->getLocalElement(globalIndex);
    if (localIndex == Teuchos::OrdinalTraits< PYTRILINOS_LOCAL_ORD >::invalid())
      return Py_BuildValue("");
    return SWIG_From_long(static_cast< long >(localIndex));
  }

  PyObject * getGlobalElement(PYTRILINOS_LOCAL_ORD localIndex)
  {
    PYTRILINOS_GLOBAL_ORD globalIndex = self->getGlobalElement(localIndex);
    if (globalIndex == Teuchos::OrdinalTraits< PYTRILINOS_GLOBAL_ORD >::invalid())
      return Py_BuildValue("");
    return SWIG_From_long(static_cast< long >(globalIndex));
  }

  PyObject * getRemoteIndexList(PyObject * globalIndexes)
  {
    // Variable initialization
    int is_new = 0;
    int type_code = PyTrilinos::NumPy_TypeCode< PYTRILINOS_GLOBAL_ORD >();
    npy_intp dims[1];
    PyObject * globalArray = NULL;
    PyObject * nodeArray   = NULL;
    PyObject * localArray  = NULL;
    PyObject * resultObj   = PyTuple_New(3);
    Teuchos::ArrayView< const PYTRILINOS_GLOBAL_ORD > globalList;
    Teuchos::ArrayView< int >                         nodeList;
    Teuchos::ArrayView< PYTRILINOS_LOCAL_ORD >        localList;
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
    type_code  = PyTrilinos::NumPy_TypeCode< PYTRILINOS_LOCAL_ORD >();
    nodeArray  = PyArray_SimpleNew(1, dims, NPY_INT);
    localArray = PyArray_SimpleNew(1, dims, type_code);
    PyTuple_SET_ITEM(resultObj, 0, nodeArray );
    PyTuple_SET_ITEM(resultObj, 1, localArray);

    // Initialize the Teuchos::ArrayViews
    globalList = Teuchos::ArrayView< const PYTRILINOS_GLOBAL_ORD >
      ((const PYTRILINOS_GLOBAL_ORD*) array_data(globalArray),
       (Teuchos::ArrayView< PYTRILINOS_GLOBAL_ORD >::size_type) dims[0]);
    nodeList = Teuchos::ArrayView< int >
      ((int*) array_data(nodeArray),
       (Teuchos::ArrayView< int >::size_type) dims[0]);
    localList = Teuchos::ArrayView< PYTRILINOS_LOCAL_ORD >
      ((PYTRILINOS_LOCAL_ORD*) array_data(localArray),
       (Teuchos::ArrayView< PYTRILINOS_LOCAL_ORD >::size_type) dims[0]);

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

%ignore Tpetra::Map< PYTRILINOS_LOCAL_ORD, PYTRILINOS_GLOBAL_ORD, DefaultNodeType >::Map;
%ignore Tpetra::Map< PYTRILINOS_LOCAL_ORD, PYTRILINOS_GLOBAL_ORD, DefaultNodeType >::getLocalElement;
%ignore Tpetra::Map< PYTRILINOS_LOCAL_ORD, PYTRILINOS_GLOBAL_ORD, DefaultNodeType >::getGlobalElement;
%ignore Tpetra::Map< PYTRILINOS_LOCAL_ORD, PYTRILINOS_GLOBAL_ORD, DefaultNodeType >::getRemoteIndexList;
%ignore Tpetra::Map< PYTRILINOS_LOCAL_ORD, PYTRILINOS_GLOBAL_ORD, DefaultNodeType >::getMyGlobalIndices;

%include "Tpetra_Map_decl.hpp"

// N.B.: Technically, the third template argument in the two SWIG
// directives below are redundant, because it is the same as the
// default template argument.  But SWIG is much more acurate when
// comparing types when all template arguments are specified.
%teuchos_rcp(Tpetra::Map< PYTRILINOS_LOCAL_ORD,
                          PYTRILINOS_GLOBAL_ORD,
                          DefaultNodeType >)
%template(Map_default) Tpetra::Map< PYTRILINOS_LOCAL_ORD,
                                    PYTRILINOS_GLOBAL_ORD,
                                    DefaultNodeType >;
%pythoncode
{
Map = Map_default
}
%inline
%{
  typedef Tpetra::Map< PYTRILINOS_LOCAL_ORD,
                       PYTRILINOS_GLOBAL_ORD,
                       DefaultNodeType > DefaultMapType;
%}

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
template <class LO = DefaultLOType,
          class GO = DefaultGOype,
          class NT = DefaultNodeType>
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
%teuchos_rcp(Tpetra::Details::Transfer< PYTRILINOS_LOCAL_ORD,
                                        PYTRILINOS_GLOBAL_ORD,
                                        DefaultNodeType >)
%template(Transfer_default) Tpetra::Details::Transfer< PYTRILINOS_LOCAL_ORD,
                                                       PYTRILINOS_GLOBAL_ORD,
                                                       DefaultNodeType >;

///////////////////////////
// Tpetra Export support //
///////////////////////////
%include "Tpetra_Export_decl.hpp"
%teuchos_rcp(Tpetra::Export< PYTRILINOS_LOCAL_ORD,
                             PYTRILINOS_GLOBAL_ORD,
                             DefaultNodeType >)
%template(Export_default) Tpetra::Export< PYTRILINOS_LOCAL_ORD,
                                          PYTRILINOS_GLOBAL_ORD,
                                          DefaultNodeType >;
%pythoncode
{
Export = Export_default
}

///////////////////////////
// Tpetra Import support //
///////////////////////////
%include "Tpetra_Import_decl.hpp"
%teuchos_rcp(Tpetra::Import< PYTRILINOS_LOCAL_ORD,
                             PYTRILINOS_GLOBAL_ORD,
                             DefaultNodeType >)
%template(Import_default) Tpetra::Import< PYTRILINOS_LOCAL_ORD,
                                          PYTRILINOS_GLOBAL_ORD,
                                          DefaultNodeType >;
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
// The refactor is making Tpetra::DistObject difficult for SWIG to
// parse, and so I provide a simplified prototype of the class here
namespace Tpetra
{
template < class Packet,
           class LocalOrdinal = DefaultLOType,
           class GlobalOrdinal = DefaultGOType,
           class Node = DefaultNodeType >
class DistObject :
    virtual public SrcDistObject,
    virtual public Teuchos::Describable
{
public:
  typedef typename Kokkos::Details::ArithTraits<Packet>::val_type packet_type;
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef DeviceType execution_space;
  typedef Kokkos::Compat::KokkosDeviceWrapperNode<execution_space> node_type;
  typedef Map<local_ordinal_type, global_ordinal_type, Node> map_type;
  explicit DistObject(const Teuchos::RCP<const DefaultMapType>& map);
  DistObject(const DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>& rhs);
  virtual ~DistObject();
  void doImport(const SrcDistObject& source,
                const Import<LocalOrdinal,GlobalOrdinal,Node>& importer,
                CombineMode CM);
  void doExport(const SrcDistObject& source,
                const Export<LocalOrdinal,GlobalOrdinal,Node>& exporter,
                CombineMode CM);
  void doImport(const SrcDistObject& source,
                const Export<LocalOrdinal,GlobalOrdinal,Node>& exporter,
                CombineMode CM);
  void doExport(const SrcDistObject& source,
                const Import<LocalOrdinal,GlobalOrdinal,Node>& importer,
                CombineMode CM);
  bool isDistributed() const;
  virtual Teuchos::RCP<const DefaultMapType> getMap() const;
  void print(std::ostream &os) const;
  virtual std::string description() const;
  virtual void
  describe(Teuchos::FancyOStream &out,
           const Teuchos::EVerbosityLevel verbLevel =
             Teuchos::Describable::verbLevel_default) const;
  virtual void
  removeEmptyProcessesInPlace(const Teuchos::RCP<const DefaultMapType>& newMap);
protected:
  virtual bool checkSizes(const SrcDistObject& source) = 0;
}; // class DistObject
} // namespace Tpetra

// %ignore Tpetra::removeEmptyProcessesInPlace;
// %include "Tpetra_DistObject_decl.hpp"
// %include "Tpetra_KokkosRefactor_DistObject_decl.hpp"

////////////////////////////////
// Tpetra MultiVector support //
////////////////////////////////
// %ignore Tpetra::MultiVector::getLocalMV;
// %ignore Tpetra::MultiVector::getLocalMVNonConst;
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
%extend Tpetra::MultiVector
{
  PyObject * _extractNumPyArray() const
  {
    if (!self->isConstantStride())
    {
      PyErr_SetString(PyExc_ValueError, "_extractNumPyArrayFromTpetraMulti"
                      "Vector: MultiVector is not constant stride");
      return NULL;
    }
    npy_intp dims[2] = { (npy_intp) self->getNumVectors(),
                         (npy_intp) self->getLocalLength() };
    const Scalar * data = self->getData(0).get();
    return PyArray_SimpleNewFromData(2,
                                     dims,
                                     PyTrilinos::NumPy_TypeCode< Scalar >(),
                                     (void*)data);
  }

  PyObject * __distarray__()
  {
    return PyTrilinos::convertToDistArray(*self);
  }
}
// The refactor is making Tpetra::MultiVector difficult for SWIG to
// parse, and so I provide a simplified prototype of the class here
namespace Tpetra
{
template< class Scalar = DefaultScalarType,
          class LocalOrdinal = DefaultLOType,
          class GlobalOrdinal = DefaultGOType,
          class Node = DefaultNodeType >
class MultiVector :
    public DistObject< Scalar,
                       LocalOrdinal,
                       GlobalOrdinal,
                       Node >
{
public:
  typedef Scalar scalar_type;
  typedef typename Kokkos::Details::ArithTraits<Scalar>::val_type impl_scalar_type;
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Kokkos::Compat::KokkosDeviceWrapperNode<Node> node_type;
  typedef typename Kokkos::Details::InnerProductSpaceTraits<impl_scalar_type>::dot_type dot_type;
  typedef typename Kokkos::Details::ArithTraits<impl_scalar_type>::mag_type mag_type;
  typedef Node execution_space;
  typedef Kokkos::DualView<impl_scalar_type**,
                           Kokkos::LayoutLeft,
                           typename execution_space::execution_space> dual_view_type;
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
  MultiVector();
  MultiVector(const Teuchos::RCP<const DefaultMapType>& map,
              const size_t numVecs,
              const bool zeroOut = true);
  MultiVector(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source);
  MultiVector(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& source,
              const Teuchos::DataAccess copyOrView);
  MultiVector(const Teuchos::RCP<const DefaultMapType>& map,
              const Teuchos::ArrayView<const Scalar>& A,
              const size_t LDA,
              const size_t NumVectors);
  MultiVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
              const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> >&ArrayOfPtrs,
              const size_t NumVectors);
  // MultiVector(const Teuchos::RCP<const DefaultMapType>& map,
  //             const dual_view_type& view);
  // MultiVector(const Teuchos::RCP<const DefaultMapType>& map,
  //             const typename dual_view_type::t_dev& d_view);
  // MultiVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
  //             const dual_view_type& view,
  //             const dual_view_type& origView);
  // MultiVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
  //             const dual_view_type& view,
  //             const Teuchos::ArrayView<const size_t>& whichVectors);
  // MultiVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
  //             const dual_view_type& view,
  //             const dual_view_type& origView,
  //             const Teuchos::ArrayView<const size_t>& whichVectors);
  template <class Node2>
  Teuchos::RCP< MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node2 > >
  clone(const Teuchos::RCP< Node2 > &node2) const;
  virtual ~MultiVector();
  void replaceGlobalValue(GlobalOrdinal globalRow,
                          size_t col,
                          const impl_scalar_type& value);
  template<typename T>
  typename std::enable_if<! std::is_same<T, impl_scalar_type>::value && std::is_convertible<T, impl_scalar_type>::value, void>::type
  replaceGlobalValue(GlobalOrdinal globalRow,
                     size_t col,
                     const T& value);
  void sumIntoGlobalValue(GlobalOrdinal globalRow,
                          size_t col,
                          const impl_scalar_type& value);
  template<typename T>
  typename std::enable_if<! std::is_same<T, impl_scalar_type>::value && std::is_convertible<T, impl_scalar_type>::value, void>::type
  sumIntoGlobalValue(GlobalOrdinal globalRow,
                     size_t col,
                     const T& value);
  void replaceLocalValue(LocalOrdinal localRow,
                         size_t col,
                         const impl_scalar_type& value);
  template<typename T>
  typename std::enable_if<! std::is_same<T, impl_scalar_type>::value && std::is_convertible<T, impl_scalar_type>::value, void>::type
  replaceLocalValue(LocalOrdinal localRow,
                    size_t col,
                    const T& value);
  void sumIntoLocalValue(LocalOrdinal localRow,
                         size_t col,
                         const impl_scalar_type& value);
  template<typename T>
  typename std::enable_if<! std::is_same<T, impl_scalar_type>::value && std::is_convertible<T, impl_scalar_type>::value, void>::type
  sumIntoLocalValue(LocalOrdinal localRow,
                    size_t col,
                    const T& value);
  void putScalar(const Scalar &value);
  void randomize();
  void replaceMap(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map);
  void reduce();
  // MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node >&
  // operator=(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node >& source);
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node > >
  subCopy(const Teuchos::Range1D &colRng) const;
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node > >
  subCopy(const Teuchos::ArrayView<const size_t> &cols) const;
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node > >
  subView(const Teuchos::Range1D &colRng) const;
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node > >
  subView(const Teuchos::ArrayView<const size_t> &cols) const;
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node > >
  subViewNonConst(const Teuchos::Range1D &colRng);
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node > >
  subViewNonConst(const Teuchos::ArrayView<const size_t> &cols);
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node > >
  offsetViewNonConst(const Teuchos::RCP<const DefaultMapType>& subMap,
                     const size_t offset);
  Teuchos::RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node > >
  getVector(const size_t j) const;
  Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node > >
  getVectorNonConst(const size_t j);
  Teuchos::ArrayRCP<const Scalar> getData(size_t j) const;
  Teuchos::ArrayRCP<Scalar> getDataNonConst(size_t j);
  void get1dCopy(const Teuchos::ArrayView<Scalar>& A,
                 const size_t LDA) const;
  void get2dCopy(const Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> >& ArrayOfPtrs) const;
  Teuchos::ArrayRCP<const Scalar> get1dView() const;
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > get2dView() const;
  Teuchos::ArrayRCP<Scalar> get1dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > get2dViewNonConst();
  // dual_view_type getDualView() const;
  template<class TargetDeviceType>
  void sync();
  template<class TargetDeviceType>
  void modify();
  template<class TargetDeviceType>
  typename Kokkos::Impl::if_c<
    Kokkos::Impl::is_same<
      typename execution_space::memory_space,
      typename TargetDeviceType::memory_space>::value,
      typename dual_view_type::t_dev,
      typename dual_view_type::t_host>::type
  getLocalView() const;
  void dot(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node >& A,
           const Teuchos::ArrayView<dot_type>& dots) const;
  template <typename T>
  typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<dot_type, T>::value), void >::type
  dot(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
      const Teuchos::ArrayView<T> &dots) const;
  template <typename T>
  typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<dot_type, T>::value), void >::type
  dot(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node >& A,
      std::vector<T>& dots) const;
  void dot(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node >& A,
           const Kokkos::View<dot_type*, execution_space>& dots) const;
  template <typename T>
  typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<dot_type, T>::value), void >::type
  dot(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node >& A,
      const Kokkos::View<T*, execution_space>& dots) const;
  void abs(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node >& A);
  void reciprocal(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A);
  void scale(const Scalar& alpha);
  void scale(Teuchos::ArrayView<const Scalar> alpha);
  void scale(const Kokkos::View<const impl_scalar_type*, execution_space> alpha);
  void scale(const Scalar& alpha,
             const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node >& A);
  void update(const Scalar& alpha,
              const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node >& A,
              const Scalar& beta);
  void update(const Scalar& alpha,
              const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node >& A,
              const Scalar& beta,
              const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node >& B,
              const Scalar& gamma);
  void norm1(const Kokkos::View<mag_type*, execution_space>& norms) const;
  template <typename T>
  typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<mag_type, T>::value), void >::type
  norm1(const Kokkos::View<T*, execution_space>& norms) const;
  void norm1(const Teuchos::ArrayView<mag_type>& norms) const;
  template <typename T>
  typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<mag_type,T>::value), void >::type
  norm1(const Teuchos::ArrayView<T>& norms) const;
  void norm2(const Kokkos::View<mag_type*, execution_space>& norms) const;
  template<typename T>
  typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<mag_type, T>::value), void >::type
  norm2(const Kokkos::View<T*, execution_space>& norms) const;
  void norm2(const Teuchos::ArrayView<mag_type>& norms) const;
  template <typename T>
  typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<mag_type,T>::value), void >::type
  norm2(const Teuchos::ArrayView<T>& norms) const;
  void normInf(const Kokkos::View<mag_type*, execution_space>& norms) const;
  template<typename T>
  typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<mag_type, T>::value), void >::type
  normInf(const Kokkos::View<T*, execution_space>& norms) const;
  void normInf(const Teuchos::ArrayView<mag_type>& norms) const;
  // typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<mag_type,T>::value), void >::type
  // normInf(const Teuchos::ArrayView<T>& norms) const;
  // void normWeighted(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node >& weights,
  //                   const Teuchos::ArrayView<mag_type>& norms) const;
  // template <typename T>
  // typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<mag_type,T>::value), void >::type
  // normWeighted(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node >& weights,
  //              const Teuchos::ArrayView<T>& norms) const;
  void meanValue(const Teuchos::ArrayView<impl_scalar_type>& means) const;
  template <typename T>
  typename Kokkos::Impl::enable_if<! Kokkos::Impl::is_same<impl_scalar_type, T>::value, void>::type
  meanValue(const Teuchos::ArrayView<T>& means) const;
  void multiply(Teuchos::ETransp transA,
                Teuchos::ETransp transB,
                const Scalar& alpha,
                const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node >& A,
                const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node >& B,
                const Scalar& beta);
  void elementWiseMultiply(Scalar scalarAB,
                           const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node >& A,
                           const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node >& B,
                           Scalar scalarThis);
  size_t getNumVectors() const;
  size_t getLocalLength() const;
  global_size_t getGlobalLength() const;
  size_t getStride() const;
  bool isConstantStride() const;
  virtual std::string description() const;
  virtual void
  describe(Teuchos::FancyOStream& out,
           const Teuchos::EVerbosityLevel verbLevel =
           Teuchos::Describable::verbLevel_default) const;
  virtual void
  removeEmptyProcessesInPlace(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& newMap);
  void setCopyOrView(const Teuchos::DataAccess copyOrView);
  Teuchos::DataAccess getCopyOrView() const;
  void assign(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node >& src);
};  // class MultiVector
}   // namespace Tpetra
// %include "Tpetra_MultiVector_decl.hpp"
// %include "Tpetra_KokkosRefactor_MultiVector_decl.hpp"
%pythoncode
%{
  def MultiVector_getattr(self, name):
      if name == "array":
          a = self._extractNumPyArray()
          self.__dict__["array"] = a
          return a
      elif name == "shape":
          return self.array.shape
      elif name == "dtype":
          return self.array.dtype
      else:
          raise AttributeError("'%s' not an attribute of MultiVector" % name)
  def MultiVector_setattr(self, name, value):
      if name in ("array", "shape", "dtype"):
          raise AttributeError("Cannot change MultiVector '%s' attribute", name)
      else:
          self.__dict__[name] = value
  def MultiVector_getitem(self,i):
      if isinstance(i,int):
          return self.getVectorNonConst(i)
      else:
          return self.array.__getitem__(i)
  def upgradeMultiVectorClass(cls):
      cls.__getattr__ = MultiVector_getattr
      cls.__setattr__ = MultiVector_setattr
      cls.__getitem__ = MultiVector_getitem
      cls.__setitem__ = lambda self, i, v: self.array.__setitem__(i,v)
      cls.__len__     = lambda self: self.array.__len__()
      cls.__str__     = lambda self: self.array.__str__()
      cls.copy        = lambda self: cls(self)
      class_array_add_math(cls)
      class_array_add_comp(cls)

%}

///////////////////////////
// Tpetra Vector support //
///////////////////////////
// The refactor is making Tpetra::Vector difficult for SWIG to
// parse, and so I provide a simplified prototype of the class here
%feature("notabstract") Tpetra::Vector;
%extend Tpetra::Vector
{
  PyObject * _extractNumPyArray() const
  {
    npy_intp dims[] = { (npy_intp) self->getLocalLength() };
    const Scalar * data = self->getData().get();
    return PyArray_SimpleNewFromData(1,
                                     dims,
                                     PyTrilinos::NumPy_TypeCode< Scalar >(),
                                     (void*)data);
  }

  PyObject * __distarray__()
  {
    return PyTrilinos::convertToDistArray(*self);
  }
}
namespace Tpetra
{
template< class Scalar = DefaultScalarType,
          class LocalOrdinal = DefaultLOType,
          class GlobalOrdinal = DefaultGOType,
          class Node = DefaultNodeType >
class Vector :
    public MultiVector< Scalar,
                        LocalOrdinal,
                        GlobalOrdinal,
                        Node >
{
public:
  typedef Scalar scalar_type;
  typedef typename base_type::impl_scalar_type impl_scalar_type;
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef typename base_type::node_type node_type;
  typedef typename base_type::dot_type dot_type;
  typedef typename base_type::mag_type mag_type;
  typedef typename base_type::dual_view_type dual_view_type;
  typedef typename base_type::map_type map_type;
  // explicit Vector(const Teuchos::RCP<const map_type>& map,
  //                 const bool zeroOut = true);
  explicit Vector(const Teuchos::RCP<const DefaultMapType >& map,
                  const bool zeroOut = true);
  Vector(const Vector<Scalar, LocalOrdinal, GlobalOrdinal,Node >& source);
  Vector(const Vector<Scalar, LocalOrdinal, GlobalOrdinal,Node >& source,
         const Teuchos::DataAccess copyOrView);

  // This constructor is giving me the following error: "Error, an
  // attempt has been made to dereference the underlying object from a
  // weak smart pointer object where the underling object has already
  // been deleted since the strong count has already gone to zero."  I
  // don't currently think it is a wrapper problem, but I could be
  // wrong.
  // Vector(const Teuchos::RCP<const DefaultMapType>& map,
  //        const Teuchos::ArrayView<const Scalar>& A);

  // I don't yet support Kokkos::DualView, so these constructors are
  // not yet appropriate.
  // Vector(const Teuchos::RCP<const DefaultMapType>& map,
  //        const dual_view_type& view);
  // Vector(const Teuchos::RCP<const DefaultMapType>& map,
  //        const dual_view_type& view,
  //        const dual_view_type& origView);

  virtual ~Vector();
  template <class Node2>
  Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node2> >
  clone(const Teuchos::RCP<Node2>& node2);
  void replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value);
  void sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value);
  void replaceLocalValue(LocalOrdinal myRow, const Scalar &value);
  void sumIntoLocalValue(LocalOrdinal myRow, const Scalar &value);
  using MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::get1dCopy;
  void get1dCopy(const Teuchos::ArrayView<Scalar>& A) const;
  using MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getDataNonConst;
  Teuchos::ArrayRCP<Scalar> getDataNonConst() { return getDataNonConst(0); }
  using MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getData;
  Teuchos::ArrayRCP<const Scalar> getData() const { return getData(0); }
  // Teuchos::RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<Node>, false> >
  // offsetView(const Teuchos::RCP<const map_type>& subMap,
  //            const size_t offset) const;
  // Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<Node>, false> >
  // offsetViewNonConst(const Teuchos::RCP<const map_type>& subMap,
  //                    const size_t offset);
  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dot;
  // dot_type dot(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<Node>, false>& y) const;
  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm1;
  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm2;
  mag_type norm2() const;
  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normInf;
  mag_type normInf() const;
  // using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normWeighted;
  // mag_type
  // normWeighted(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<Node>, false>& weights) const;
  using MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, false>::meanValue;
  Scalar meanValue() const;
  virtual std::string description() const;
  virtual void
  describe(Teuchos::FancyOStream& out,
           const Teuchos::EVerbosityLevel verbLevel =
           Teuchos::Describable::verbLevel_default) const;
};  // class Vector
}   // namespace Tpetra
// %ignore Tpetra::Vector::getLocalMV;
// %ignore Tpetra::Vector::getLocalMVNonConst;
// %warnfilter(302) Tpetra::createVectorFromView;
// %include "Tpetra_Vector_decl.hpp"
// %include "Tpetra_KokkosRefactor_Vector_decl.hpp"
%pythoncode
%{
  def Vector_getattr(self, name):
      if name == "array":
          a = self._extractNumPyArray()
          self.__dict__["array"] = a
          return a
      elif name == "shape":
          return self.array.shape
      elif name == "dtype":
          return self.array.dtype
      else:
          raise AttributeError("'%s' not an attribute of Vector" % name)
  def Vector_setattr(self, name, value):
      if name in ("array", "shape", "dtype"):
          raise AttributeError("Cannot change Vector '%s' attribute", name)
      else:
          self.__dict__[name] = value
  def upgradeVectorClass(cls):
      cls.__getattr__ = Vector_getattr
      cls.__setattr__ = Vector_setattr
      cls.__getitem__ = lambda self, i: self.array.__getitem__(i)
      cls.__setitem__ = lambda self, i, v: self.array.__setitem__(i,v)
      cls.__len__     = lambda self: self.array.__len__()
      cls.__str__     = lambda self: self.array.__str__()
      cls.copy        = lambda self: cls(self)
      class_array_add_math(cls)
      class_array_add_comp(cls)

%}

/////////////////////////////////////
// Explicit template instantiation //
/////////////////////////////////////
//
// Macro names:
//     CLASS        Class name (DistObject, Vector, ...)
//     SCALAR       C/C++ scalar name (int, long long, ...)
//     SCALAR_NAME  Suffix name (int, longlong, ...)
//
%define %tpetra_class(CLASS, SCALAR, SCALAR_NAME)
    %warnfilter(315) Tpetra::CLASS< SCALAR,
                                    PYTRILINOS_LOCAL_ORD,
                                    PYTRILINOS_GLOBAL_ORD,
                                    DefaultNodeType >;
    %template(CLASS ## _ ## SCALAR_NAME) Tpetra::CLASS< SCALAR,
                                                        PYTRILINOS_LOCAL_ORD,
                                                        PYTRILINOS_GLOBAL_ORD,
                                                        DefaultNodeType >;
%enddef

%define %tpetra_scalars(SCALAR, SCALAR_NAME)
    %teuchos_rcp(Tpetra::DistObject< SCALAR,
                                     PYTRILINOS_LOCAL_ORD,
                                     PYTRILINOS_GLOBAL_ORD,
                                     DefaultNodeType >)
    %tpetra_class(DistObject, SCALAR, SCALAR_NAME)

    %teuchos_rcp_dap(PyTrilinos::convertPythonToTpetraMultiVector< SCALAR >,
                     Tpetra::MultiVector< SCALAR,
                                          PYTRILINOS_LOCAL_ORD,
                                          PYTRILINOS_GLOBAL_ORD,
                                          DefaultNodeType >)
    %tpetra_class(MultiVector, SCALAR, SCALAR_NAME)
    %pythoncode
    %{
      upgradeMultiVectorClass(MultiVector_ ## SCALAR_NAME)
    %}

    %teuchos_rcp_dap(PyTrilinos::convertPythonToTpetraVector< SCALAR >,
                     Tpetra::Vector< SCALAR,
                                     PYTRILINOS_LOCAL_ORD,
                                     PYTRILINOS_GLOBAL_ORD,
                                     DefaultNodeType >)
    %tpetra_class(Vector, SCALAR, SCALAR_NAME)
    %pythoncode
    %{
      upgradeVectorClass(Vector_ ## SCALAR_NAME)
    %}
%enddef

//////////////////////////////////////////////
// Concrete scalar types for Tpetra classes //
//////////////////////////////////////////////
%tpetra_scalars(int       , int   )
%tpetra_scalars(long long , long  )
%tpetra_scalars(double    , double)

/////////////////////////////////////////////////////
// Python code that consolidates templated classes //
/////////////////////////////////////////////////////
%pythoncode
%{
  def MultiVector(*args, **kwargs):
    dtype = None
    if len(args) > 0:
      try:
        dtype = str(args[0].dtype)
      except AttributeError:
        pass
    dtype = kwargs.get("dtype", dtype)
    if dtype is None: dtype = "int64"
    if type(dtype) == str:
      dtype = numpy.dtype(dtype)
    if dtype.type is numpy.int32:
      result = MultiVector_int(*args)
    elif dtype.type is numpy.int64:
      result = MultiVector_long(*args)
    #elif dtype.type is numpy.float32:
    #  result = MultiVector_float(*args)
    elif dtype.type is numpy.float64:
      result = MultiVector_double(*args)
    else:
      raise TypeError("Unsupported or unrecognized dtype = %s" %
                      str(dtype))
    return result

  def Vector(*args, **kwargs):
    dtype = None
    if len(args) > 0:
      try:
        dtype = str(args[0].dtype)
      except AttributeError:
        pass
    dtype = kwargs.get("dtype", dtype)
    if dtype is None: dtype = "int64"
    if type(dtype) == str:
      dtype = numpy.dtype(dtype)
    if dtype.type is numpy.int32:
      result = Vector_int(*args)
    elif dtype.type is numpy.int64:
      result = Vector_long(*args)
    #elif dtype.type is numpy.float32:
    #  result = Vector_float(*args)
    elif dtype.type is numpy.float64:
      result = Vector_double(*args)
    else:
      raise TypeError("Unsupported or unrecognized dtype = %s" %
                      str(dtype))
    return result
%}

// Turn off exception handling
%exception;
