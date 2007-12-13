// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2005) Sandia Corporation
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

%{
#include "Teuchos_VerbosityLevel.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_LabeledObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ReductionOp.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#include "Teuchos_CommHelpers.hpp"
%}

/////////////////////////////////////
// Teuchos::VerbosityLevel support //
/////////////////////////////////////
%rename(verbosityLevelToString) Teuchos::toString;
%include "Teuchos_VerbosityLevel.hpp"

///////////////////////////////////
// Teuchos::FancyOStream support //
///////////////////////////////////
%include "Teuchos_FancyOStream.hpp"

////////////////////////////////////
// Teuchos::LabeledObject support //
////////////////////////////////////
%include "Teuchos_LabeledObject.hpp"

//////////////////////////////////
// Teuchos::Describable support //
//////////////////////////////////
%include "Teuchos_Describable.hpp"

//////////////////////////////////
// Teuchos::ReductionOp support //
//////////////////////////////////
%include "Teuchos_ReductionOp.hpp"

///////////////////////////
// Teuchos::Comm support //
///////////////////////////
%extend Teuchos::Comm
{
  PyObject * broadcast(int rootRank, PyObject * buffer) const
  {
    PyArrayObject * bufferArray = obj_to_array_no_conversion(buffer, NPY_NOTYPE);
    if (!bufferArray || !require_contiguous(bufferArray)) return NULL;
    Ordinal bytes = static_cast<Ordinal>(PyArray_NBYTES(bufferArray));
    char * bufferVals = (char*) array_data(bufferArray);
    self->broadcast(rootRank, bytes, bufferVals);
    return Py_BuildValue("");
  }
  PyObject * gatherAll(PyObject * sendBuffer) const
  {
    int     is_new_object  = 0;
    Ordinal sendBytes      = 0;
    Ordinal recvBytes      = 0;
    int     recvNd         = 0;
    int     type           = 0;
    char *  sendBufferVals = NULL;
    char *  recvBufferVals = NULL;
    PyArrayObject * recvBufferArray = NULL;
    PyArrayObject * sendBufferArray =
      obj_to_array_contiguous_allow_conversion(sendBuffer, NPY_NOTYPE, &is_new_object);
    if (!sendBufferArray) goto fail;
    sendBytes = static_cast<Ordinal>(PyArray_NBYTES(sendBufferArray));
    recvBytes = sendBytes * self->getSize();
    recvNd    = array_numdims(sendBufferArray) + 1;
    type      = array_type(sendBufferArray);
    { // Scope this to make recvDims temporary
      npy_intp * recvDims = new npy_intp[recvNd];
      recvDims[0] = self->getSize();
      for (int i=1; i<recvNd; ++i) recvDims[i] = array_size(sendBufferArray,i-1);
      recvBufferArray = (PyArrayObject*) PyArray_SimpleNew(recvNd, recvDims, type);
      if (!recvBufferArray) goto fail;
      delete [] recvDims;
    }
    sendBufferVals = (char*) array_data(sendBufferArray);
    recvBufferVals = (char*) array_data(recvBufferArray);
    self->gatherAll(sendBytes, sendBufferVals, recvBytes, recvBufferVals);
    if (is_new_object && sendBufferArray) Py_DECREF(sendBufferArray);
    return PyArray_Return(recvBufferArray);
  fail:
    if (is_new_object && sendBufferArray) Py_DECREF(sendBufferArray);
    return NULL;
  }
  PyObject * reduceAll(Teuchos::EReductionType reductOp, PyObject * sendBuffer) const
  {
    int     is_new_object = 0;
    Ordinal bytes         = 0;
    int     type          = 0;
    char * sendBufferVals    = NULL;
    char * globalReductsVals = NULL;
    PyArrayObject * globalReductsArray = NULL;
    PyArrayObject * sendBufferArray    =
      obj_to_array_contiguous_allow_conversion(sendBuffer, NPY_NOTYPE, &is_new_object);
    if (!sendBufferArray) goto fail;
    bytes = static_cast<Ordinal>(PyArray_NBYTES(sendBufferArray));
    type  = array_type(sendBufferArray);
    globalReductsArray = (PyArrayObject*)
      PyArray_SimpleNew(array_numdims(sendBufferArray),
			array_dimensions(sendBufferArray), type);
    PyArray_FILLWBYTE(globalReductsArray, 0);
    sendBufferVals    = (char*) array_data(sendBufferArray);
    globalReductsVals = (char*) array_data(globalReductsArray);
    Teuchos::reduceAll(*self, reductOp, bytes, sendBufferVals, globalReductsVals);
    if (is_new_object && sendBufferArray) Py_DECREF(sendBufferArray);
    return PyArray_Return(globalReductsArray);
  fail:
    if (is_new_object && sendBufferArray) Py_DECREF(sendBufferArray);
    return NULL;
  }
  PyObject * scan(Teuchos::EReductionType reductOp, PyObject * sendBuffer) const
  {
    int     is_new_object = 0;
    Ordinal bytes         = 0;
    int     type          = 0;
    char * sendBufferVals  = NULL;
    char * scanReductsVals = NULL;
    PyArrayObject * scanReductsArray = NULL;
    PyArrayObject * sendBufferArray  =
      obj_to_array_contiguous_allow_conversion(sendBuffer, NPY_NOTYPE, &is_new_object);
    if (!sendBufferArray) goto fail;
    bytes = static_cast<Ordinal>(PyArray_NBYTES(sendBufferArray));
    type  = array_type(sendBufferArray);
    scanReductsArray = (PyArrayObject*)
      PyArray_SimpleNew(array_numdims(sendBufferArray),
			array_dimensions(sendBufferArray), type);
    PyArray_FILLWBYTE(scanReductsArray, 0);
    sendBufferVals    = (char*) array_data(sendBufferArray);
    scanReductsVals = (char*) array_data(scanReductsArray);
    Teuchos::scan(*self, reductOp, bytes, sendBufferVals, scanReductsVals);
    if (is_new_object && sendBufferArray) Py_DECREF(sendBufferArray);
    return PyArray_Return(scanReductsArray);
  fail:
    if (is_new_object && sendBufferArray) Py_DECREF(sendBufferArray);
    return NULL;
  }
  std::string __str__() const
  {
    return self->description();
  }
}
%ignore Teuchos::Comm::broadcast;
%ignore Teuchos::Comm::gatherAll;
%ignore Teuchos::Comm::reduceAll;
%ignore Teuchos::Comm::scan;
%include "Teuchos_Comm.hpp"
%template(Comm_long) Teuchos::Comm<long>;
%pythoncode %{
Comm = Comm_long
%}

/////////////////////////////////
// Teuchos::SerialComm support //
/////////////////////////////////
%ignore Teuchos::SerialComm::broadcast;
%ignore Teuchos::SerialComm::gatherAll;
%ignore Teuchos::SerialComm::reduceAll;
%ignore Teuchos::SerialComm::scan;
%include "Teuchos_DefaultSerialComm.hpp"
%template(SerialComm_long) Teuchos::SerialComm<long>;
%pythoncode %{
SerialComm = SerialComm_long
PyComm = SerialComm
%}

//////////////////////////////////
// Teuchos::Comm Helper support //
//////////////////////////////////
%rename(reductionTypeToString) Teuchos::toString;
%include "Teuchos_CommHelpers.hpp"
%template(rank_long   ) Teuchos::rank<long>;
%template(size_long   ) Teuchos::size<long>;
%template(barrier_long) Teuchos::barrier<long>;
%pythoncode %{
rank    = rank_long
size    = size_long
barrier = barrier_long

def broadcast(comm, rootRank, buffer):
  """
  broadcast(Teuchos.Comm comm, int rootRank, numpy.ndarray buffer)

  Broadcast the contents of buffer from processor rootRank to all of the other
  processors.
  """
  comm.broadcast(rootRank, buffer)

def gatherAll(comm, buffer):
  """
  gatherAll(Teuchos.Comm comm, buffer) -> numpy.ndarray

  Gather the contents of buffer to all of the processors.
  """
  return comm.gatherAll(buffer)

def reduceAll(comm, reductOp, buffer):
  """
  Reduce the contents of buffer according to the operation designated by
  reductOp on all of the processors.
  """
  return comm.reduceAll(reductOp, buffer)

def scan(comm, reductOp, buffer):
  """
  Return the scan of the contents of buffer according to the operation
  designated by reductOp on each of the processors.
  """
  return comm.scan(reductOp, buffer)
%}

