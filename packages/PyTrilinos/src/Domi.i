// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//     Domi: Multi-dimensional Distributed Linear Algebra Services
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

%define %domi_docstring
"
PyTrilinos.Domi is the python interface to the Trilinos structured,
multi-dimensional, distrivuted linear algebra servise package Domi:

    http://trilinos.sandia.gov/packages/domi

Domi supports the structured decomposition of structured vectors
(arrays), maps, and communicators.  It also supports the conversion of
these objects to Epetra and Tpetra Vectors, MultiVectors, and Maps
(including as views where possible), so that they can be used with
other Trilinos solver technologies.
"
%enddef

%module(package   = "PyTrilinos",
        autodoc   = "1",
        docstring = %domi_docstring) Domi

%{
// System include files
#include <iostream>
#include <sstream>
#include <vector>

// Teuchos include files
#include "PyTrilinos_Teuchos_Headers.hpp"

// Epetra include files
#ifdef HAVE_EPETRA
#include "PyTrilinos_Epetra_Headers.hpp"
#endif

// Tpetra include files
#ifdef HAVE_TPETRA
#include "PyTrilinos_Tpetra_Headers.hpp"
#endif

// Domi include files
#include "PyTrilinos_Domi_Headers.hpp"

%}

// Auto-documentation feature
%feature("autodoc", "1");

%include "PyTrilinos_config.h"

// Domi enumerated types support
#undef  PACKAGE_BUGREPORT
%ignore PACKAGE_BUGREPORT;
#undef  PACKAGE_NAME
%ignore PACKAGE_NAME;
#undef  PACKAGE_STRING
%ignore PACKAGE_STRING;
#undef  PACKAGE_TARNAME
%ignore PACKAGE_TARNAME;
#undef  PACKAGE_VERSION
%ignore PACKAGE_VERSION;
%include "Domi_config.h"
%include "Domi_ConfigDefs.hpp"

// Include Domi documentation
%include "Domi_dox.i"

// SWIG library include files
%include "stl.i"

// Include the NumPy typemaps
%include "numpy.i"
%pythoncode
{
import numpy
}

// Include the standard exception handlers
%include "exception.i"
%include "Domi_exceptions.i"

// External Teuchos interface imports
%import "Teuchos.i"
%include "Teuchos_Array.i"
%teuchos_array_typemaps(Domi::dim_type , NPY_INT )
%teuchos_array_typemaps(Domi::size_type, NPY_LONG)

// External Epetra interface imports
#ifdef HAVE_EPETRA
%import "Epetra.i"
#endif

// External Tpetra interface imports
#ifdef HAVE_TPETRA
%import "Tpetra.i"

// Define shortcuts for the default Tpetra template types
%inline
%{
  typedef Tpetra::Details::DefaultTypes::scalar_type         DefaultScalarType;
  typedef Tpetra::Details::DefaultTypes::local_ordinal_type  DefaultLOType;
  typedef Tpetra::Details::DefaultTypes::global_ordinal_type DefaultGOType;
  typedef Tpetra::Details::DefaultTypes::node_type           DefaultNodeType;

%}
#endif

// General exception handling
%feature("director:except")
{
  if ($error != NULL)
  {
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
  catch(int errCode)
  {
    PyErr_Format(PyExc_RuntimeError, "Error code = %d\nSee stderr for details",
                 errCode);
    SWIG_fail;
  }
  catch(Swig::DirectorException &e)
  {
    SWIG_fail;
  }
  SWIG_CATCH_DOMIEXCEPT
  SWIG_CATCH_STDEXCEPT
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

// General ignore directives
%ignore operator<<;
%ignore Domi::operator<<;
%ignore operator==;
%ignore Domi::operator==;
%ignore operator!=;
%ignore Domi::operator!=;

////////////////////////////
// Domi Utilities support //
////////////////////////////
%ignore Domi::remove_const;
%ignore Domi::computeStrides;
%ignore Domi::computeSize;
%ignore Domi::regularizeCommDims;
%ignore Domi::computeCommIndexes;
%ignore Domi::computePeriodic;
%ignore Domi::splitStringOfIntsWithCommas;
#ifdef HAVE_MPI
%ignore Domi::mpiType;
%ignore Domi::mpiOrder;
#endif
%include "Domi_Utils.hpp"

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

//////////////////////////
// Domi Version support //
//////////////////////////
%include "Domi_Version.hpp"
%pythoncode
%{
  __version__ = Version().split()[2]
%}

/////////////////////////////////////
// Domi getValidParameters support //
/////////////////////////////////////
%include "Domi_getValidParameters.hpp"

////////////////////////
// Domi Slice support //
////////////////////////
// Note that we do not wrap the Domi::Slice class, but rather provide
// typemaps for conversion to and from Python slice objects.
%include "Domi_Slice.i"
// We do, however, need to %import Domi_Slice.hpp, so that SWIG knows
// that Slice is in the Domi namespace.
%import "Domi_Slice.hpp"

//////////////////////////
// Domi MDArray support //
//////////////////////////
%ignore Domi::swap;
%include "Domi_MDArray.i"

// /////////////////////////
// // Domi MDComm support //
// /////////////////////////
%teuchos_rcp(Domi::MDComm)
%ignore Domi::MDComm::operator=;
%include "Domi_MDComm.hpp"
%extend Domi::MDComm
{
  Domi::MDComm __getitem__(PyObject * indexes)
  {
    // If 'indexes' is not a sequence, it might be an integer or
    // slice.  So wrap it in a tuple, and we'll check its type below.
    if (!PySequence_Check(indexes))
    {
      PyObject * newIndexes = Py_BuildValue("(N)", indexes);
      indexes = newIndexes;
    }

    // Get the number of indexes in the sequence.  If this is larger
    // than the number of dimensions of the MDComm, then cap it at
    // that value.
    Py_ssize_t numIndexes = PySequence_Size(indexes);
    if (numIndexes > self->numDims()) numIndexes = self->numDims();

    // Initialize the new MDComm as a copy of this MDComm
    Domi::MDComm newMdComm(*self);

    // 'domiAxis' will be the index for the new MDComm as we construct
    // it.  'axis' will be the index for the sequence of indexes.
    // These can diverge as the new MDComm is constructed.
    int domiAxis = 0;
    for (Py_ssize_t axis = 0; axis < numIndexes; ++axis)
    {
      if (newMdComm.onSubcommunicator())
      {
        PyObject * index = PySequence_GetItem(indexes, axis);
        if (PyInt_Check(index))
        {
          int axisRank = (int) PyInt_AsLong(index);
          newMdComm = Domi::MDComm(newMdComm, domiAxis, axisRank);
          // Do not increment domiAxis, because the new MDComm has one
          // fewer dimension!
        }
        else if (PySlice_Check(index))
        {
          PySliceObject * pySlice = (PySliceObject*) index;
          Py_ssize_t commDim = (Py_ssize_t) newMdComm.getCommDim(domiAxis);
          Domi::Slice slice = PyTrilinos::convertToDomiSlice(pySlice, commDim);
          newMdComm = Domi::MDComm(newMdComm, domiAxis, slice);
          domiAxis++;
        }
        else
        {
          PyErr_SetString(PyExc_TypeError, "Argument type error for "
                          "Domi.MDComm __getitem__.  Argument must be a "
                          "sequence of integers and/or slices");
          throw PyTrilinos::PythonException();
        }
        Py_DECREF(index);
      }
    }
    return newMdComm;
  }
}

////////////////////////
// Domi MDMap support //
////////////////////////
%extend Domi::MDMap< Node >
{
  Domi::MDMap< Node > __getitem__(PyObject * indexes)
  {
    // If 'indexes' is not a sequence, it might be an integer or
    // slice.  So wrap it in a tuple, and we'll check its type below.
    if (!PySequence_Check(indexes))
    {
      PyObject * newIndexes = Py_BuildValue("(N)", indexes);
      indexes = newIndexes;
    }

    // Get the number of indexes in the sequence.  If this is larger
    // than the number of dimensions of the MDMap, then cap it at
    // that value.
    Py_ssize_t numIndexes = PySequence_Size(indexes);
    if (numIndexes > self->numDims()) numIndexes = self->numDims();

    // Initialize the new MDMap as a copy of this MDMap
    Domi::MDMap< Node > newMdMap(*self);

    // 'domiAxis' will be the index for the new MDMap as we construct
    // it.  'axis' will be the index for the sequence of indexes.
    // These can diverge as the new MDMap is constructed.
    int domiAxis = 0;
    for (Py_ssize_t axis = 0; axis < numIndexes; ++axis)
    {
      if (newMdMap.onSubcommunicator())
      {
        PyObject * index = PySequence_GetItem(indexes, axis);
        if (PyInt_Check(index))
        {
          int axisRank = (int) PyInt_AsLong(index);
          newMdMap = Domi::MDMap< Node >(newMdMap, domiAxis, axisRank);
          // Do not increment domiAxis, because the new MDMap has one
          // fewer dimension!
        }
        else if (PySlice_Check(index))
        {
          PySliceObject * pySlice = (PySliceObject*) index;
          Py_ssize_t dim = (Py_ssize_t) newMdMap.getGlobalDim(domiAxis);
          Domi::Slice slice = PyTrilinos::convertToDomiSlice(pySlice, dim);
          newMdMap = Domi::MDMap< Node >(newMdMap, domiAxis, slice);
          domiAxis++;
        }
        else
        {
          PyErr_SetString(PyExc_TypeError, "Argument type error for "
                          "Domi.MDMap __getitem__.  Argument must be a "
                          "sequence of integers and/or slices");
          throw PyTrilinos::PythonException();
        }
        Py_DECREF(index);
      }
    }
    return newMdMap;
  }
}
%include "Domi_MDMap.hpp"
#ifdef HAVE_TPETRA
%template(getTpetraMap) Domi::MDMap::getTpetraMap< PYTRILINOS_LOCAL_ORD,
                                                   PYTRILINOS_GLOBAL_ORD,
                                                   DefaultNodeType >;
%template(getTpetraAxisMap) Domi::MDMap::getTpetraAxisMap< PYTRILINOS_LOCAL_ORD,
                                                           PYTRILINOS_GLOBAL_ORD,
                                                           DefaultNodeType >;
#endif
%teuchos_rcp(Domi::MDMap< Domi::DefaultNode::DefaultNodeType >)
%template(MDMap_default) Domi::MDMap< Domi::DefaultNode::DefaultNodeType >;
%pythoncode
{
MDMap = MDMap_default
}

///////////////////////////
// Domi MDVector support //
///////////////////////////
%extend Domi::MDVector
{
  Domi::MDVector< Scalar, Node > __getitem__(PyObject * indexes)
  {
    // If 'indexes' is not a sequence, it might be an integer or
    // slice.  So wrap it in a tuple, and we will check its type below.
    if (!PySequence_Check(indexes))
    {
      PyObject * newIndexes = Py_BuildValue("(N)", indexes);
      indexes = newIndexes;
    }

    // Get the number of indexes in the sequence.  If this is larger
    // than the number of dimensions of the MDVector, then cap it at
    // that value.
    Py_ssize_t numIndexes = PySequence_Size(indexes);
    if (numIndexes > self->numDims()) numIndexes = self->numDims();

    // Initialize the new MDVector as a view of this MDVector
    Domi::MDVector< Scalar, Node > newMdVector(*self, Teuchos::View);

    // 'domiAxis' will be the index for the new MDVector as we construct
    // it.  'axis' will be the index for the sequence of indexes.
    // These can diverge as the new MDVector is constructed.
    int domiAxis = 0;
    for (Py_ssize_t axis = 0; axis < numIndexes; ++axis)
    {
      if (newMdVector.onSubcommunicator())
      {
        PyObject * index = PySequence_GetItem(indexes, axis);
        if (PyInt_Check(index))
        {
          int axisRank = (int) PyInt_AsLong(index);
          newMdVector = Domi::MDVector< Scalar, Node >(newMdVector,
                                                       domiAxis,
                                                       axisRank);
          // Do not increment domiAxis, because the new MDVector has one
          // fewer dimension!
        }
        else if (PySlice_Check(index))
        {
          PySliceObject * pySlice = (PySliceObject*) index;
          Py_ssize_t dim = (Py_ssize_t) newMdVector.getGlobalDim(domiAxis);
          Domi::Slice slice = PyTrilinos::convertToDomiSlice(pySlice, dim);
          newMdVector = Domi::MDVector< Scalar, Node >(newMdVector,
                                                       domiAxis,
                                                       slice);
          domiAxis++;
        }
        else
        {
          PyErr_SetString(PyExc_TypeError, "Argument type error for "
                          "Domi.MDVector __getitem__.  Argument must be a "
                          "sequence of integers and/or slices");
          throw PyTrilinos::PythonException();
        }
        Py_DECREF(index);
      }
    }
    return newMdVector;
  }

  Domi::MDArrayView< Scalar > getData(bool includePadding = true)
  {
    return self->getDataNonConst(includePadding);
  }

  PyObject * __distarray__()
  {
    return PyTrilinos::convertToDistArray(*self);
  }

  // PyObject * dtype()
  // {
  //   return PyArray_TypeObjectFromType(PyTrilinos::NumPy_TypeCode< Scalar >());
  // }
}
%ignore Domi::MDVector::operator=;
%ignore Domi::MDVector::operator[];
%ignore Domi::MDVector::getDataNonConst(bool includePadding = true);
%ignore Domi::MDVector::getData(bool includePadding = true) const;
%include "Domi_MDVector.hpp"
%teuchos_rcp(Domi::MDVector< int      , Domi::DefaultNode::DefaultNodeType >)
%teuchos_rcp(Domi::MDVector< long long, Domi::DefaultNode::DefaultNodeType >)
%teuchos_rcp(Domi::MDVector< double   , Domi::DefaultNode::DefaultNodeType >)
//%teuchos_rcp(Domi::MDVector< float    , Domi::DefaultNode::DefaultNodeType >)
%pythoncode
%{
  def MDVector_getattr(self, name):
      if name == "array":
          a = self.getData()
          self.__dict__["array"] = a
          return a
      elif name == "shape":
          return self.array.shape
      elif name == "dtype":
          return self.array.dtype
      else:
          raise AttributeError("'%s' not an attribute of MDVector" % name)
  def MDVector_setattr(self, name, value):
      if name in ("array", "shape", "dtype"):
          raise AttributeError("Cannot change MDVector '%s' attribute", name)
      else:
          self.__dict__[name] = value
  def upgradeMDVectorClass(cls):
      cls.__getattr__ = MDVector_getattr
      cls.__setattr__ = MDVector_setattr
      cls.__setitem__ = lambda self, i, v: self.array.__setitem__(i,v)
      cls.__len__     = lambda self: self.array.__len__()
      cls.__str__     = lambda self: self.array.__str__()
      cls.copy        = lambda self: cls(self)
      class_array_add_math(cls)
      class_array_add_comp(cls)

%}
#ifdef HAVE_TPETRA
%template(getTpetraVectorView)
    Domi::MDVector::getTpetraVectorView< PYTRILINOS_LOCAL_ORD,
                                         PYTRILINOS_GLOBAL_ORD,
                                         DefaultNodeType >;
%template(getTpetraVectorCopy)
    Domi::MDVector::getTpetraVectorCopy< PYTRILINOS_LOCAL_ORD,
                                         PYTRILINOS_GLOBAL_ORD,
                                         DefaultNodeType >;
%template(getTpetraMultiVectorView)
    Domi::MDVector::getTpetraMultiVectorView< PYTRILINOS_LOCAL_ORD,
                                              PYTRILINOS_GLOBAL_ORD,
                                              DefaultNodeType >;
%template(getTpetraMultiVectorCopy)
    Domi::MDVector::getTpetraMultiVectorCopy< PYTRILINOS_LOCAL_ORD,
                                              PYTRILINOS_GLOBAL_ORD,
                                              DefaultNodeType >;
#endif
%template(MDVector_int   )
  Domi::MDVector< int      , Domi::DefaultNode::DefaultNodeType >;
%pythoncode
%{
  upgradeMDVectorClass(MDVector_int)
%}
%template(MDVector_long  )
  Domi::MDVector< long long, Domi::DefaultNode::DefaultNodeType >;
%pythoncode
%{
  upgradeMDVectorClass(MDVector_long)
%}
%template(MDVector_double)
  Domi::MDVector< double   , Domi::DefaultNode::DefaultNodeType >;
%pythoncode
%{
  upgradeMDVectorClass(MDVector_double)
%}
// %template(MDVector_float )
//   Domi::MDVector< float    , Domi::DefaultNode::DefaultNodeType >;
// %pythoncode
// %{
//   upgradeMDVectorClass(MDVector_float)
// %}

////////////////////////////
// from_DistArray support //
////////////////////////////
%inline
{
template< class Scalar >
Teuchos::RCP< Domi::MDVector< Scalar, Domi::DefaultNode::DefaultNodeType > >
from_DistArray(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
               PyObject * distArrayObj)
{
  // if (!PyObject_HasAttrString(distArrayObj, "__distarray__"))
  // {
  //   PyErr_SetString(PyExc_ValueError, "Object does not have '__distarray__'"
  //                   " method");
  //   throw PyTrilinos::PythonException();
  // }
  // PyObject * distarray = PyObject_GetAttrString(distArrayObj, "__distarray__");
  // PyTrilinos::DistArrayProtocol dap(distarray);
  PyTrilinos::DistArrayProtocol dap(distArrayObj);
  return PyTrilinos::convertToMDVector< Scalar >(teuchosComm, dap);
}
}
%template(from_DistArray_int   ) from_DistArray< int       >;
%template(from_DistArray_long  ) from_DistArray< long long >;
%template(from_DistArray_double) from_DistArray< double    >;
//%template(from_DistArray_float ) from_DistArray< float     >;
%pythoncode
%{
def from_DistArray(comm, distarray):
    protocol = distarray.__distarray__()
    dtype = protocol["buffer"].dtype
    if dtype.type is numpy.int32:
        return from_DistArray_int(comm, protocol)
    elif dtype.type is numpy.int64:
        return from_DistArray_long(comm, protocol)
    #elif dtype.type is numpy.float32:
    #    return from_DistArray_float(comm, protocol)
    elif dtype.type is numpy.float64:
        return from_DistArray_double(comm, protocol)
    else:
        raise TypeError("Unsupported or unrecognized dtype = %s" % str(dtype))
%}

%pythoncode
%{
# My first cut at MDVector was to make it a wrapper class.  I now prefer the
# idea of MDVector being a factory function.  There is a bit of non-trivial work
# that went into the wrapper class, though, so I comment it out rather than
# delete it.  Who knows, I may want to go back to it some day...

# class MDVector(object):
#     def __init__(self, *args, **kwargs):
#         dtype       = kwargs.get("dtype"      , "int64")
#         zeroOut     = kwargs.get("zeroOut"    , False  )
#         leadingDim  = kwargs.get("leadingDim" , 0      )
#         trailingDim = kwargs.get("trailingDim", 0      )
#         if type(dtype) == str:
#             dtype = numpy.dtype(dtype)
# 
#         # Factory for arg is MDMap
#         if isinstance(args[0], MDMap):
#             if dtype.type is numpy.int32:
#                 self._vector = MDVector_int(args[0],
#                                             leadingDim,
#                                             trailingDim,
#                                             zeroOut)
#             elif dtype.type is numpy.int64:
#                 self._vector = MDVector_long(args[0],
#                                              leadingDim,
#                                              trailingDim,
#                                              zeroOut)
#             elif dtype.type is numpy.float32:
#                 self._vector = MDVector_float(args[0],
#                                               leadingDim,
#                                               trailingDim,
#                                               zeroOut)
#             elif dtype.type is numpy.float64:
#                 self._vector = MDVector_double(args[0],
#                                                leadingDim,
#                                                trailingDim,
#                                                zeroOut)
#             else:
#                 raise TypeError("Unsupported or unrecognized dtype = %s" %
#                                 str(dtype))
# 
#         # Factory for arg is DistArray
#         elif hasattr(arg, '__distarray__'):
#             self._vector = from_DistArray(*args)
# 
#         self.dtype = dtype
# 
#     def __getattribute__(self, name):
#         if name in ('__class__', '__dir__', '__getitem__', '_vector', 'dtype'):
#             return object.__getattribute__(self, name)
#         return getattr(object.__getattribute__(self, '_vector'), name)
# 
#     def __dir__(self):
#         return sorted(set(dir(self._vector) + dir(MDVector)))
# 
#     def __getitem__(self, args):
#         return self._vector.__getitem__(args)

def MDVector(*args, **kwargs):
    dtype = None
    if len(args) > 0:
        try:
            dtype = str(args[0].dtype())
            if dtype == "int": dtype = "i"
        except AttributeError:
            pass
    dtype = kwargs.get("dtype", dtype)
    if dtype is None: dtype = "int64"
    if type(dtype) == str:
        dtype = numpy.dtype(dtype)
    if dtype.type is numpy.int32:
        result = MDVector_int(*args)
    elif dtype.type is numpy.int64:
        result = MDVector_long(*args)
    #elif dtype.type is numpy.float32:
    #    result = MDVector_float(*args)
    elif dtype.type is numpy.float64:
        result = MDVector_double(*args)
    else:
        raise TypeError("Unsupported or unrecognized dtype = %s" %
                        str(dtype))
    return result

%}

// Turn off the exception handling
%exception;
