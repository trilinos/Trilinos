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

// Local includes
#include "PyTrilinos_Tpetra_Util.hpp"
#include "PyTrilinos_PythonException.hpp"

namespace PyTrilinos
{

////////////////////////////////////////////////////////////////////////

PyObject *
convertToDimData(const Teuchos::RCP< const Tpetra::Map< PYTRILINOS_LOCAL_ORD,
                                                        PYTRILINOS_GLOBAL_ORD > > & tm,
                 int   extraDim)
{
  // Initialization
  PyObject * dim_data    = NULL;
  PyObject * dim_dict    = NULL;
  PyObject * indices     = NULL;
  npy_intp   dims        = 1;
  int        currentDim  = 0;
  Teuchos::ArrayView< const PYTRILINOS_GLOBAL_ORD > nodeBlockIDs;

  // Get the Teuchos::Comm
  Teuchos::RCP< const Teuchos::Comm< int > > comm = tm->getComm();

  // Get the number of dimensions.  The vector data constitutes one
  // dimension.  If the extraDim is greater than one, then that
  // constitutes a second possible dimension.
  Py_ssize_t ndim = 1;
  if (extraDim > 1) ndim += 1;

  // Allocate the dim_data return object, which is a tuple of length
  // ndim
  dim_data = PyTuple_New(ndim);
  if (!dim_data) goto fail;

  // If we have an extra dimension argument greater than one, then
  // define a dimension associated with the multiple vectors
  if (extraDim > 1)
  {
    dim_dict = PyDict_New();
    if (!dim_dict) goto fail;
    if (PyDict_SetItemString(dim_dict,
                             "dist_type",
                             PyString_FromString("n")) == -1) goto fail;

    if (PyDict_SetItemString(dim_dict,
                             "size",
                             PyInt_FromLong(extraDim)) == -1) goto fail;
    // This function steals a reference from dim_dict
    PyTuple_SET_ITEM(dim_data, currentDim, dim_dict);
    currentDim += 1;
  }

  // Allocate the dimension data dictionary for the distributed points
  // and assign the common key values
  dim_dict = PyDict_New();
  if (!dim_dict) goto fail;
  if (PyDict_SetItemString(dim_dict,
                           "size",
                           PyInt_FromLong(tm->getGlobalNumElements())) == -1)
    goto fail;
  if (PyDict_SetItemString(dim_dict,
                           "proc_grid_size",
                           PyInt_FromLong(comm->getSize())) == -1) goto fail;
  if (PyDict_SetItemString(dim_dict,
                           "proc_grid_rank",
                           PyInt_FromLong(comm->getRank())) == -1) goto fail;

  // Determine the type of the dimension data, either block "b" or
  // unstructured "u", set the "dist_type" key and other keys required
  // according to dist_type.
  if (tm->isContiguous())
  {
    // isContiguous() == true corresponds to DistArray Protocol
    // dist_type == "b" (block)
    if (PyDict_SetItemString(dim_dict,
                             "dist_type",
                             PyString_FromString("b")) == -1) goto fail;
    if (PyDict_SetItemString(dim_dict,
                             "start",
                             PyInt_FromLong(tm->getMinGlobalIndex())) == -1)
      goto fail;
    if (PyDict_SetItemString(dim_dict,
                             "stop",
                             PyInt_FromLong(tm->getMaxGlobalIndex()+1)) == -1)
      goto fail;
  }
  else
  {
    // isContiguous() == false corresponds to DistArray Protocol
    // dist_type == "u" (unstructured)
    if (PyDict_SetItemString(dim_dict,
                             "dist_type",
                             PyString_FromString("u")) == -1) goto fail;
    dims    = tm->getNodeElementList().size();
    nodeBlockIDs = tm->getNodeElementList();
    indices = PyArray_SimpleNewFromData(1,
                                        &dims,
                                        NPY_LONG,
                                        (void*)nodeBlockIDs.getRawPtr());
    if (!indices) goto fail;
    if (PyDict_SetItemString(dim_dict,
                             "indices",
                             indices) == -1) goto fail;
    Py_DECREF(indices);
    indices = NULL;
  }

  // Put the dimension dictionary into the dim_data tuple
  PyTuple_SET_ITEM(dim_data, currentDim, dim_dict);

  // Return the dim_data tuple
  return dim_data;

  // Error handling
  fail:
  Py_XDECREF(dim_data);
  Py_XDECREF(dim_dict);
  Py_XDECREF(indices);
  return NULL;
}

////////////////////////////////////////////////////////////////////////

}  // Namespace PyTrilinos
