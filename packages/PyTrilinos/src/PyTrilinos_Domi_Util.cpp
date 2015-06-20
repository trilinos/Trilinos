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

#include "PyTrilinos_Domi_Util.hpp"
#include "PyTrilinos_PythonException.hpp"

////////////////////////////////////////////////////////////////////////

namespace PyTrilinos
{

////////////////////////////////////////////////////////////////////////

PyObject *
convertToPySlice(const Domi::Slice & domiSlice)
{
  PyObject * pyStart = NULL;
  PyObject * pyStop  = NULL;
  PyObject * pyStep  = NULL;
  PyObject * pySlice;
  if (domiSlice.start() != Domi::Slice::Default)
    pyStart = PyInt_FromLong((long)domiSlice.start());
  if (domiSlice.stop()  != Domi::Slice::Default)
    pyStop = PyInt_FromLong((long)domiSlice.stop());
  if (domiSlice.step()  != Domi::Slice::Default)
    pyStep = PyInt_FromLong((long)domiSlice.step());
  return PySlice_New(pyStart, pyStop, pyStep);
}

////////////////////////////////////////////////////////////////////////

Domi::Slice
convertToDomiSlice(PySliceObject * pySlice,
                   Py_ssize_t length)
{
  Py_ssize_t pyStart  = 0;
  Py_ssize_t pyStop   = 0;
  Py_ssize_t pyStep   = 0;
  Py_ssize_t sliceLen = 0;
  int rcode = PySlice_GetIndicesEx(pySlice, length , &pyStart ,
                                   &pyStop, &pyStep, &sliceLen);
  Domi::Slice domiSlice((Domi::dim_type) pyStart,
                        (Domi::dim_type) pyStop ,
                        (Domi::dim_type) pyStep );
  return domiSlice;
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< const Domi::MDComm >
convertToMDComm(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
                const DistArrayProtocol & distarray)
{
  // Get the number of dimensions and allocate the MDComm arrays
  int numDims = distarray.num_dims();
  Teuchos::Array< int > commDims(numDims);
  Teuchos::Array< int > periodic(numDims);

  // Fill in the MDComm arrays
  for (int axis = 0; axis < numDims; ++axis)
  {
    commDims[axis] = distarray.proc_grid_size(axis);
    periodic[axis] = distarray.periodic(axis) ? 1 : 0;
  }

  // Return the result
  return Teuchos::rcp(new Domi::MDComm(teuchosComm, commDims, periodic));
}

////////////////////////////////////////////////////////////////////////

Teuchos::RCP< const Domi::MDMap<> >
convertToMDMap(const Teuchos::RCP< const Teuchos::Comm< int > > teuchosComm,
               const DistArrayProtocol & distarray)
{
  // Get the equivalent MDComm.
  Teuchos::RCP< const Domi::MDComm > mdComm =
    convertToMDComm(teuchosComm, distarray);

  // Get the number of dimensions and check the types of the
  // distributions
  int numDims = mdComm->numDims();
  for (int axis = 0; axis < numDims; ++axis)
  {
    DistributionType distType = distarray.dist_type(axis);

    // Check the distribution type
    if (distType == NONE)
    {
      PyErr_Format(PyExc_ValueError,
                   "'dist_type' for axis %d = 'NONE' is invalid",
                   axis);
      throw PythonException();
    }
    if (distType == CYCLIC)
    {
      PyErr_Format(PyExc_ValueError,
                   "'dist_type' for axis %d = 'CYCLIC' is invalid",
                   axis);
      throw PythonException();
    }
    if (distType == UNSTRUCTURED)
    {
      PyErr_Format(PyExc_ValueError,
                   "'dist_type' for axis %d = 'UNSTRUCTURED' is invalid",
                   axis);
      throw PythonException();
    }
  }

  // Allocate the MDMap constructor arrays
  Teuchos::Array< Domi::Slice > myGlobalBounds(numDims);
  Teuchos::Array< PaddingType > padding(numDims);

  // Fill in the MDMap constructor arrays
  for (int axis = 0; axis < numDims; ++axis)
  {
    int start = distarray.start(axis);
    int stop  = distarray.stop( axis);
    myGlobalBounds[axis] = Domi::ConcreteSlice(start, stop);
    padding[axis]  = distarray.padding( axis);
  }

  // Get the buffer layout
  Domi::Layout layout =
    PyArray_IS_C_CONTIGUOUS((PyArrayObject*) distarray.buffer()) ?
      Domi::C_ORDER : Domi::FORTRAN_ORDER;

  // Return the result
  return Teuchos::rcp(new Domi::MDMap<>(mdComm,
                                        myGlobalBounds,
                                        padding,
                                        layout));
}

////////////////////////////////////////////////////////////////////////

PyObject * convertToDimData(const Teuchos::RCP< const Domi::MDMap<> > mdMap)
{
  Py_ssize_t numDims = mdMap->numDims();
  PyObject * dimData = PyTuple_New(numDims);
  Teuchos::RCP< const Teuchos::Comm< int > > comm = mdMap->getTeuchosComm();
  Domi::Slice bounds;
  int loPad;
  int hiPad;
  PyObject * periodic;
  PyObject * dimDict;
  for (int axis = 0; axis < numDims; ++axis)
  {
    bounds   = mdMap->getGlobalRankBounds(axis,false);
    loPad    = mdMap->getLowerPadSize(axis);
    hiPad    = mdMap->getUpperPadSize(axis);
    periodic = mdMap->isPeriodic(axis) ? Py_True : Py_False;
    dimDict  = PyDict_New();
    if (PyDict_SetItemString(dimDict, "dist_type",
                             Py_BuildValue("s", "b")) == -1)
      goto fail;
    if (PyDict_SetItemString(dimDict, "size",
                             Py_BuildValue("i", mdMap->getLocalDim(axis,true)))
        == -1) goto fail;
    if (PyDict_SetItemString(dimDict, "proc_grid_size",
                             Py_BuildValue("i", mdMap->getCommDim(axis)))
        == -1) goto fail;
    if (PyDict_SetItemString(dimDict, "proc_grid_rank",
                             Py_BuildValue("i", mdMap->getCommIndex(axis)))
        == -1) goto fail;
    if (PyDict_SetItemString(dimDict, "start",
                             Py_BuildValue("i", bounds.start()-loPad)) == -1)
      goto fail;
    if (PyDict_SetItemString(dimDict, "stop",
                             Py_BuildValue("i", bounds.stop()+hiPad)) == -1)
      goto fail;
    if (PyDict_SetItemString(dimDict, "padding",
                             Py_BuildValue("(ii)", loPad, hiPad)) == -1)
      goto fail;
    Py_INCREF(periodic);
    if (PyDict_SetItemString(dimDict, "periodic", periodic) == -1)
      goto fail;
    PyTuple_SET_ITEM(dimData, axis, dimDict);
  }
  return dimData;

  fail:
  if (PyDict_Check(dimDict)) PyDict_Clear(dimDict);
  for (int axis = 0; axis < numDims; ++axis)
  {
    dimDict = PyTuple_GET_ITEM(dimData, axis);
    if (dimDict) PyDict_Clear(dimDict);
  }
  return NULL;
}

////////////////////////////////////////////////////////////////////////

}
