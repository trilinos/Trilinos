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

// Include the Python prototypes
#include "Python.h"

// Include the NumPy and Python headers
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"

// PyTrilinos include
#include "PyTrilinos_config.h"
#include "PyTrilinos_NumPy_Util.hpp"

// Teuchos includes
#ifdef HAVE_TEUCHOS
#include "Teuchos_RCP.hpp"
#endif

// Epetra includes
//#include "Tpetra_BlockMap.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"

////////////////////////////////////////////////////////////////////////

namespace PyTrilinos
{

////////////////////////////////////////////////////////////

// Given a Tpetra::Map, return a Python dimension data object, which
// is a tuple of Python dimension data dictionaries that describe the
// Tpetra::Map, consistent with the DistArray Protocol.  The extraDim
// argument is to allow for the multiple vectors of a
// Tpetra::MultiVector.  If an error occurs, return NULL.  Note that a
// Tpetra::Map with variable element sizes is currently not supported
// and results in an error.
PyObject *
convertTpetraMapToDimData(const Tpetra::Map< long, long > & tm,
                          int   extraDim=1);

////////////////////////////////////////////////////////////

// Given a Tpetra::BlockMap, return a Python dimension data object,
// which is a tuple of Python dimension data dictionaries that
// describe the Tpetra::BlockMap, consistent with the DistArray
// Protocol.  The extraDim argument is to allow for the multiple
// vectors of a Tpetra::MultiVector.  If an error occurs, return NULL.
// Note that a Tpetra::BlockMap with variable element sizes is
// currently not supported and results in an error.
// PyObject *
// convertTpetraBlockMapToDimData(const Tpetra::BlockMap< long, long > & tbm,
//                                int   extraDim=1);

////////////////////////////////////////////////////////////

// Given a Tpetra::MultiVector, return a Python dictionary consistent
// with the DistArray Protocol.  If an error occurs, return NULL.
template< class Scalar >
PyObject *
convertTpetraMultiVectorToDAP(const Tpetra::MultiVector< Scalar > & tmv);

////////////////////////////////////////////////////////////////////////
// *** Implementations ***
////////////////////////////////////////////////////////////////////////

template< class Scalar >
PyObject *
convertTpetraMultiVectorToDAP(const Tpetra::MultiVector< Scalar > & tmv)
{
  // Initialization
  PyObject   * dap      = NULL;
  PyObject   * dim_data = NULL;
  PyObject   * dim_dict = NULL;
  PyObject   * size     = NULL;
  PyObject   * buffer   = NULL;
  Py_ssize_t   ndim     = 1;
  npy_intp     dims[3];
  Teuchos::ArrayRCP< Scalar > data;

  // Get the underlying Tpetra::Map< long, long >
  Teuchos::RCP< const Tpetra::Map< long, long > > tm = tmv.getMap();

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
  dim_data = convertTpetraBlockMapToDimData(*tm,
                                            tmv.getNumVectors());
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
    size = PyDict_GetItemString(dim_dict, "size");
    if (!size) goto fail;
    dims[i] = PyInt_AsLong(size);
    if (PyErr_Occurred()) goto fail;
  }
  data = tmv.getDataNonConst();
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
