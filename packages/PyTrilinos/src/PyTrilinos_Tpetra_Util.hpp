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

#ifndef PYTRILINOS_TPETRA_UTIL_HPP
#define PYTRILINOS_TPETRA_UTIL_HPP

// Standard includes
#include <string.h>

// Include the Python prototypes
#include "Python.h"

// Include the NumPy and Python headers
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"

// PyTrilinos includes
#include "PyTrilinos_config.h"
#include "Python3Compat.hpp"
#include "PyTrilinos_NumPy_Util.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Tpetra includes
#include "Tpetra_Map.hpp"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "Tpetra_MultiVector.hpp"
#pragma GCC diagnostic pop

////////////////////////////////////////////////////////////////////////

namespace PyTrilinos
{

// Given a Tpetra::Map, return a Python dimension data object, which
// is a tuple of Python dimension data dictionaries that describe the
// Tpetra::Map, consistent with the DistArray Protocol.  The extraDim
// argument is to allow for the multiple vectors of a
// Tpetra::MultiVector.  If an error occurs, return NULL.  Note that a
// Tpetra::Map with variable element sizes is currently not supported
// and results in an error.
PyObject *
convertToDimData(const Teuchos::RCP< const Tpetra::Map< PYTRILINOS_LOCAL_ORD,
                                                        PYTRILINOS_GLOBAL_ORD > > & tm,
                 int   extraDim=1);

////////////////////////////////////////////////////////////

// Given a Tpetra::MultiVector, return a Python dictionary consistent
// with the DistArray Protocol.  If an error occurs, return NULL.
template< class Scalar >
PyObject *
convertToDistArray(const Tpetra::MultiVector< Scalar,
                                              PYTRILINOS_LOCAL_ORD,
                                              PYTRILINOS_GLOBAL_ORD > & tmv);

////////////////////////////////////////////////////////////////////////
// *** Implementations ***
////////////////////////////////////////////////////////////////////////

template< class Scalar >
PyObject *
convertToDistArray(const Tpetra::MultiVector< Scalar,
                                              PYTRILINOS_LOCAL_ORD,
                                              PYTRILINOS_GLOBAL_ORD > & tmv)
{
  // Initialization
  PyObject   * dap       = NULL;
  PyObject   * dim_data  = NULL;
  PyObject   * dim_dict  = NULL;
  PyObject   * dist_type = NULL;
  PyObject   * start     = NULL;
  PyObject   * stop      = NULL;
  PyObject   * indices   = NULL;
  PyObject   * buffer    = NULL;
  Py_ssize_t   ndim      = 1;
  npy_intp     dims[3];
  Teuchos::ArrayRCP< const Scalar > data;

  // Get the underlying Tpetra::Map< PYTRILINOS_LOCAL_ORD, PYTRILINOS_GLOBAL_ORD >
  Teuchos::RCP< const Tpetra::Map< PYTRILINOS_LOCAL_ORD,
                                   PYTRILINOS_GLOBAL_ORD > > tm = tmv.getMap();

  // Allocate the DistArray Protocol object and set the version key
  // value
  dap = PyDict_New();
  if (!dap) goto fail;
  if (PyDict_SetItemString(dap,
                           "__version__",
                           PyString_FromString("0.9.0")) == -1) goto fail;

  // Get the Dimension Data and the number of dimensions.  If the
  // underlying Tpetra::BlockMap has variable element sizes, an error
  // will be detected here.
  dim_data = convertToDimData(tm, tmv.getNumVectors());
  if (!dim_data) goto fail;
  ndim = PyTuple_Size(dim_data);

  // Assign the Dimension Data key value.
  if (PyDict_SetItemString(dap,
                           "dim_data",
                           dim_data) == -1) goto fail;

  // Extract the buffer dimensions from the Dimension Data, construct
  // the buffer and assign the buffer key value
  for (Py_ssize_t i = 0; i < ndim; ++i)
  {
    dim_dict = PyTuple_GetItem(dim_data, i);
    if (!dim_dict) goto fail;
    dist_type = PyDict_GetItemString(dim_dict, "dist_type");
    if (!dist_type) goto fail;
    if (strcmp(convertPyStringToChar(dist_type), "b") == 0)
    {
      start = PyDict_GetItemString(dim_dict, "start");
      if (!start) goto fail;
      stop = PyDict_GetItemString(dim_dict, "stop");
      if (!stop) goto fail;
      dims[i] = PyInt_AsLong(stop) - PyInt_AsLong(start);
      if (PyErr_Occurred()) goto fail;
    }
    else if (strcmp(convertPyStringToChar(dist_type), "u") == 0)
    {
      indices = PyDict_GetItemString(dim_dict, "indices");
      if (!indices) goto fail;
      dims[i] = PyArray_DIM((PyArrayObject*)indices,0);
      if (PyErr_Occurred()) goto fail;
    }
    else
    {
      PyErr_Format(PyExc_ValueError,
                   "Unsupported distribution type '%s'",
                   convertPyStringToChar(dist_type));
      goto fail;
    }
  }
  data = tmv.getData(0);
  buffer = PyArray_SimpleNewFromData(ndim,
                                     dims,
                                     NumPy_TypeCode< Scalar >(),
                                     (void*)data.getRawPtr());
  if (!buffer) goto fail;
  if (PyDict_SetItemString(dap,
                           "buffer",
                           buffer) == -1) goto fail;

  // Return the DistArray Protocol object
  return dap;

  // Error handling
  fail:
  Py_XDECREF(dap);
  Py_XDECREF(dim_data);
  Py_XDECREF(dim_dict);
  Py_XDECREF(buffer);
  return NULL;
}

}

#endif

#if defined(PyTrilinos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The PyTrilinos package is deprecated"
#endif
#endif

