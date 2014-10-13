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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%define %domi_docstring
"
PyTrilinos.Domi is the python interface to the Trilinos structured,
multi-dimensional, distrivuted linear algebra servise package Domi:

    http://trilinos.sandia.gov/preCopyrightTrilinos/domi

Domi supports the structured decomposition of structured vectors
(arrays), maps, and communicators.  It also supports the conversion of
these objects to Epetra and Tpetra Vectors, MultiVectors, and Maps
(including as views where possible), so that they can be used with
other Trilinos solver technologies.
"
%enddef

%module(package = "PyTrilinos",
        autodoc   = "1",
	docstring = %domi_docstring) Domi

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

// Configuration includes
#include "PyTrilinos_config.h"
#ifdef HAVE_SYS_TIME_H
#undef HAVE_SYS_TIME_H
#endif
#ifdef HAVE_INTTYPES_H
#undef HAVE_INTTYPES_H
#endif
#ifdef HAVE_STDINT_H
#undef HAVE_STDINT_H
#endif
#include "Domi_ConfigDefs.hpp"

// Teuchos includes
#include "Teuchos_CommHelpers.hpp"

#ifdef HAVE_EPETRA
// Epetra includes
#include "Epetra_LAPACK.h"
#include "Epetra_Time.h"
#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_SerialDistributor.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_OffsetIndex.h"
#include "Epetra_LocalMap.h"
#include "Epetra_IntVector.h"
#include "Epetra_FEVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_InvOperator.h"
#include "Epetra_MapColoring.h"
#include "Epetra_SerialDenseSVD.h"
#include "Epetra_SerialDenseSolver.h"
#endif

// Domi includes
#include "Domi_Utils.hpp"
#include "Domi_Version.hpp"
#include "Domi_getValidParameters.hpp"
#include "Domi_Slice.hpp"
#include "Domi_MDArrayView.hpp"
#include "Domi_MDArray.hpp"
#include "Domi_MDArrayRCP.hpp"
#include "Domi_MDComm.hpp"
#include "Domi_MDMap.hpp"
#include "Domi_MDVector.hpp"

// PyTrilinos includes
#include "PyTrilinos_Teuchos_Util.hpp"
#include "Epetra_NumPyIntSerialDenseVector.hpp"
#include "Epetra_NumPyIntSerialDenseMatrix.hpp"
#include "Epetra_NumPySerialDenseVector.hpp"
#include "Epetra_NumPySerialDenseMatrix.hpp"
#include "Epetra_NumPySerialSymDenseMatrix.hpp"
#include "Epetra_NumPyIntVector.hpp"
#include "Epetra_NumPyVector.hpp"
#include "Epetra_NumPyFEVector.hpp"
#include "Epetra_NumPyMultiVector.hpp"
#include "PyTrilinos_Domi_Util.hpp"

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
//%include "Domi_dox.i"

// SWIG library includes
%include "stl.i"

// Include the NumPy typemaps
%include "numpy.i"
%pythoncode
{
import numpy
}

// Include the standard exception handlers
%include "exception.i"

// External Teuchos interface imports
%import "Teuchos.i"
%include "Teuchos_Array.i"
%teuchos_array_typemaps(Domi::dim_type , NPY_INT )
%teuchos_array_typemaps(Domi::size_type, NPY_LONG)

// External Epetra interface imports
#ifdef HAVE_EPETRA
%import "Epetra.i"
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
%ignore Domi::MDMap::getTpetraMap;
%ignore Domi::MDMap::getTpetraAxisMap;
%extend Domi::MDMap
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
%teuchos_rcp(Domi::MDMap< Kokkos::DefaultNode::DefaultNodeType >)
%template(MDMap_default) Domi::MDMap< Kokkos::DefaultNode::DefaultNodeType >;
%pythoncode
{
MDMap = MDMap_default
}

///////////////////////////
// Domi MDVector support //
///////////////////////////
%ignore Domi::MDVector::operator=;
%ignore Domi::MDVector::operator[];
%ignore Domi::MDVector::getTpetraVectorView;
%ignore Domi::MDVector::getTpetraMultiVectorView;
%ignore Domi::MDVector::getTpetraVectorCopy;
%ignore Domi::MDVector::getTpetraMultiVectorCopy;
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
}
%ignore Domi::MDVector::getDataNonConst;
%ignore Domi::MDVector::getData;

%include "Domi_MDVector.hpp"
%teuchos_rcp(Domi::MDVector< int    >)
%template(MDVector_int   ) Domi::MDVector< int    >;
%teuchos_rcp(Domi::MDVector< long   >)
%template(MDVector_long  ) Domi::MDVector< long   >;
%teuchos_rcp(Domi::MDVector< float  >)
%template(MDVector_float ) Domi::MDVector< float  >;
%teuchos_rcp(Domi::MDVector< double >)
%template(MDVector_double) Domi::MDVector< double >;

////////////////////////////
// from_DistArray support //
////////////////////////////
%inline
{
template< class Scalar >
Teuchos::RCP< Domi::MDVector< Scalar > >
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
%template(from_DistArray_int   ) from_DistArray< int    >;
%template(from_DistArray_long  ) from_DistArray< long   >;
%template(from_DistArray_float ) from_DistArray< float  >;
%template(from_DistArray_double) from_DistArray< double >;
%pythoncode
{
def from_DistArray(comm, distarray):
    protocol = distarray.__distarray__()
    dtype = protocol["buffer"].dtype
    if dtype.type is numpy.int32:
        return from_DistArray_int(comm, protocol)
    elif dtype.type is numpy.int64:
        return from_DistArray_long(comm, protocol)
    elif dtype.type is numpy.float32:
        return from_DistArray_float(comm, protocol)
    elif dtype.type is numpy.float64:
        return from_DistArray_double(comm, protocol)
    else:
        raise TypeError("Unsupported or unrecognized dtype = %s" % str(dtype))
}

%pythoncode
{
class MDVector(object):
    def __init__(self, *args, **kwargs):
        dtype       = kwargs.get("dtype"      , "int64")
        zeroOut     = kwargs.get("zeroOut"    , False  )
        leadingDim  = kwargs.get("leadingDim" , 0      )
        trailingDim = kwargs.get("trailingDim", 0      )
        if type(dtype) == str:
            dtype = numpy.dtype(dtype)

        # Factory for arg is MDMap
        if isinstance(args[0], MDMap):
            if dtype.type is numpy.int32:
                self._vector = MDVector_int(args[0],
                                            leadingDim,
                                            trailingDim,
                                            zeroOut)
            elif dtype.type is numpy.int64:
                self._vector = MDVector_long(args[0],
                                             leadingDim,
                                             trailingDim,
                                             zeroOut)
            elif dtype.type is numpy.float32:
                self._vector = MDVector_float(args[0],
                                              leadingDim,
                                              trailingDim,
                                              zeroOut)
            elif dtype.type is numpy.float64:
                self._vector = MDVector_double(args[0],
                                               leadingDim,
                                               trailingDim,
                                               zeroOut)
            else:
                raise TypeError("Unsupported or unrecognized dtype = %s" %
                            str(dtype))

        # Factory for arg is DistArray
        elif hasattr(arg, '__distarray__'):
            self._vector = from_DistArray(*args)

        self.__dtype = dtype

    def __getattribute__(self, name):
        if name in ('__class__', '__dir__', '__getitem__', '_vector'):
            return object.__getattribute__(self, name)
        return getattr(object.__getattribute__(self, '_vector'), name)

    def __dir__(self):
        return sorted(set(dir(self._vector) + dir(MDVector)))

    # __getitem__ has to be pulled out specifically, probably because it comes
    # from %extend, although I do not understand why
    def __getitem__(self, args):
        return self._vector.__getitem__(args)
}

// Turn off the exception handling
%exception;
