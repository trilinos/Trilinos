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
%feature("director") Teuchos::LabeledObject;
%include "Teuchos_LabeledObject.hpp"

//////////////////////////////////
// Teuchos::Describable support //
//////////////////////////////////
%feature("director") Teuchos::Describable;
%include "Teuchos_Describable.hpp"

//////////////////////////////////
// Teuchos::ReductionOp support //
//////////////////////////////////
%include "Teuchos_ReductionOp.hpp"

///////////////////////////
// Teuchos::Comm support //
///////////////////////////
%feature("autodoc",
"broadcast(self, int rootRank, numpy.ndarray buffer)

Broadcast the contents of buffer from processor rootRank to all of the
other processors.  Argument buffer must be a numpy array, so that the
broadcast can be performed in-place.  Its scalar data type can be any
numerical type supported by numpy.")
Teuchos::Comm::broadcast;
%feature("autodoc",
"gatherAll(self, buffer) -> numpy.ndarray

Gather the contents of buffer to all of the processors.  Argument
buffer can be a numpy array or any sequence that can be converted to a
numpy array.  Its scalar data type can be any numerical type supported
by numpy.  The return argument is a numpy array of the same type.")
Teuchos::Comm::gatherAll;
%feature("autodoc",
"reduceAll(EReductionType reductOp, buffer) -> numpy.ndarray

Reduce the contents of buffer according to the operation designated by
reductOp on all of the processors.  Argument reductOp can be
Teuchos.REDUCE_SUM, Teuchos.REDUCE_MAX, or Teuchos.REDUCE_MIN.
Argument buffer can be a numpy array or any sequence that can be
converted to a numpy array.  Its scalar data type can be any numerical
type supported by numpy.  The return argument is a numpy array of the
same type.")
Teuchos::Comm::reduceAll;
%feature("autodoc",
"scan(EReductionType reductOp, buffer) -> numpy.ndarray

Return the scan of the contents of buffer according to the operation
designated by reductOp on each of the processors.  Argument reductOp
can be Teuchos.REDUCE_SUM, Teuchos.REDUCE_MAX, or Teuchos.REDUCE_MIN.
Argument buffer can be a numpy array or any sequence that can be
converted to a numpy array.  Its scalar data type can be any numerical
type supported by numpy.  The return argument is a numpy array of the
same type.")
Teuchos::Comm::reduceAll;
%extend Teuchos::Comm
{
  PyObject * broadcast(int rootRank, PyObject * bcastObj) const
  {
    PyArrayObject * bcastArray = obj_to_array_no_conversion(bcastObj, NPY_NOTYPE);
    if (!bcastArray || !require_contiguous(bcastArray)) return NULL;
    Ordinal bytes = static_cast<Ordinal>(PyArray_NBYTES(bcastArray));
    char * bcastBuffer = (char*) array_data(bcastArray);
    self->broadcast(rootRank, bytes, bcastBuffer);
    return Py_BuildValue("");
  }

  PyObject * gatherAll(PyObject * sendObj) const
  {
    int     is_new_object = 0;
    Ordinal sendBytes     = 0;
    Ordinal recvBytes     = 0;
    int     recvNd        = 0;
    int     type          = 0;
    char *  sendBuffer    = NULL;
    char *  recvBuffer    = NULL;
    PyArrayObject * recvArray = NULL;
    PyArrayObject * sendArray =
      obj_to_array_contiguous_allow_conversion(sendObj, NPY_NOTYPE, &is_new_object);
    if (!sendArray) goto fail;
    sendBytes = static_cast<Ordinal>(PyArray_NBYTES(sendArray));
    recvBytes = sendBytes * self->getSize();
    recvNd    = array_numdims(sendArray) + 1;
    type      = array_type(sendArray);
    { // Scope this to make recvDims temporary
      npy_intp * recvDims = new npy_intp[recvNd];
      recvDims[0] = self->getSize();
      for (int i=1; i<recvNd; ++i) recvDims[i] = array_size(sendArray,i-1);
      recvArray = (PyArrayObject*) PyArray_SimpleNew(recvNd, recvDims, type);
      if (!recvArray) goto fail;
      delete [] recvDims;
    }
    sendBuffer = (char*) array_data(sendArray);
    recvBuffer = (char*) array_data(recvArray);
    self->gatherAll(sendBytes, sendBuffer, recvBytes, recvBuffer);
    if (is_new_object && sendArray) { Py_DECREF(sendArray); }
    return PyArray_Return(recvArray);
  fail:
    if (is_new_object && sendArray) { Py_DECREF(sendArray); }
    return NULL;
  }

  PyObject * reduceAll(Teuchos::EReductionType reductOp, PyObject * sendObj) const
  {
    int     is_new_object = 0;
    Ordinal count         = 0;
    int     type          = 0;
    void *  sendBuffer    = NULL;
    void *  globalReducts = NULL;
    PyArrayObject * globalArray = NULL;
    PyArrayObject * sendArray   =
      obj_to_array_contiguous_allow_conversion(sendObj, NPY_NOTYPE, &is_new_object);
    if (!sendArray) goto fail;
    count = static_cast<Ordinal>(PyArray_SIZE(sendArray));
    type  = array_type(sendArray);
    globalArray = (PyArrayObject*)
      PyArray_SimpleNew(array_numdims(sendArray),
			array_dimensions(sendArray), type);
    PyArray_FILLWBYTE(globalArray, 0);
    sendBuffer    = array_data(sendArray);
    globalReducts = array_data(globalArray);
    switch (type)
    {
    case NPY_BYTE:
      Teuchos::reduceAll(*self, reductOp, count, (char*)sendBuffer,
			 (char*)globalReducts); break;
    case NPY_UBYTE:
      Teuchos::reduceAll(*self, reductOp, count, (unsigned char*)sendBuffer,
			 (unsigned char*)globalReducts); break;
    case NPY_SHORT:
      Teuchos::reduceAll(*self, reductOp, count, (short*)sendBuffer,
			 (short*)globalReducts); break;
    case NPY_USHORT:
      Teuchos::reduceAll(*self, reductOp, count, (unsigned short*)sendBuffer,
			 (unsigned short*)globalReducts); break;
    case NPY_INT:
      Teuchos::reduceAll(*self, reductOp, count, (int*)sendBuffer,
			 (int*)globalReducts); break;
    case NPY_UINT:
      Teuchos::reduceAll(*self, reductOp, count, (unsigned int*)sendBuffer,
			 (unsigned int*)globalReducts); break;
    case NPY_LONG:
      Teuchos::reduceAll(*self, reductOp, count, (long*)sendBuffer,
			 (long*)globalReducts); break;
    case NPY_ULONG:
      Teuchos::reduceAll(*self, reductOp, count, (unsigned long*)sendBuffer,
			 (unsigned long*)globalReducts); break;
    case NPY_LONGLONG:
      Teuchos::reduceAll(*self, reductOp, count, (long long*)sendBuffer,
			 (long long*)globalReducts); break;
    case NPY_ULONGLONG:
      Teuchos::reduceAll(*self, reductOp, count, (unsigned long long*)sendBuffer,
			 (unsigned long long*)globalReducts); break;
    case NPY_FLOAT:
      Teuchos::reduceAll(*self, reductOp, count, (float*)sendBuffer,
			 (float*)globalReducts); break;
    case NPY_DOUBLE:
      Teuchos::reduceAll(*self, reductOp, count, (double*)sendBuffer,
			 (double*)globalReducts); break;
    default:
      PyErr_SetString(PyExc_TypeError, "reduceAll() for unsupported NumPy type");
      goto fail;
    }
    if (is_new_object && sendArray) { Py_DECREF(sendArray); }
    return PyArray_Return(globalArray);
  fail:
    if (is_new_object && sendArray) { Py_DECREF(sendArray); }
    Py_XDECREF(globalArray);
    return NULL;
  }

  PyObject * scan(Teuchos::EReductionType reductOp, PyObject * sendObj) const
  {
    int     is_new_object = 0;
    Ordinal count         = 0;
    int     type          = 0;
    void *  sendBuffer    = NULL;
    void *  scanReducts   = NULL;
    PyArrayObject * scanArray = NULL;
    PyArrayObject * sendArray =
      obj_to_array_contiguous_allow_conversion(sendObj, NPY_NOTYPE, &is_new_object);
    if (!sendArray) goto fail;
    count = static_cast<Ordinal>(PyArray_SIZE(sendArray));
    type  = array_type(sendArray);
    scanArray = (PyArrayObject*)
      PyArray_SimpleNew(array_numdims(sendArray),
			array_dimensions(sendArray), type);
    PyArray_FILLWBYTE(scanArray, 0);
    sendBuffer  = array_data(sendArray);
    scanReducts = array_data(scanArray);
    switch (type)
    {
    case NPY_BYTE:
      Teuchos::scan(*self, reductOp, count, (char*)sendBuffer,
		    (char*)scanReducts); break;
    case NPY_UBYTE:
      Teuchos::scan(*self, reductOp, count, (unsigned char*)sendBuffer,
		    (unsigned char*)scanReducts); break;
    case NPY_SHORT:
      Teuchos::scan(*self, reductOp, count, (short*)sendBuffer,
		    (short*)scanReducts); break;
    case NPY_USHORT:
      Teuchos::scan(*self, reductOp, count, (unsigned short*)sendBuffer,
		    (unsigned short*)scanReducts); break;
    case NPY_INT:
      Teuchos::scan(*self, reductOp, count, (int*)sendBuffer,
		    (int*)scanReducts); break;
    case NPY_UINT:
      Teuchos::scan(*self, reductOp, count, (unsigned int*)sendBuffer,
		    (unsigned int*)scanReducts); break;
    case NPY_LONG:
      Teuchos::scan(*self, reductOp, count, (long*)sendBuffer,
		    (long*)scanReducts); break;
    case NPY_ULONG:
      Teuchos::scan(*self, reductOp, count, (unsigned long*)sendBuffer,
		    (unsigned long*)scanReducts); break;
    case NPY_LONGLONG:
      Teuchos::scan(*self, reductOp, count, (long long*)sendBuffer,
		    (long long*)scanReducts); break;
    case NPY_ULONGLONG:
      Teuchos::scan(*self, reductOp, count, (unsigned long long*)sendBuffer,
		    (unsigned long long*)scanReducts); break;
    case NPY_FLOAT:
      Teuchos::scan(*self, reductOp, count, (float*)sendBuffer,
		    (float*)scanReducts); break;
    case NPY_DOUBLE:
      Teuchos::scan(*self, reductOp, count, (double*)sendBuffer,
		    (double*)scanReducts); break;
    default:
      PyErr_SetString(PyExc_TypeError, "scan() for unsupported NumPy type");
      goto fail;
    }
    if (is_new_object && sendArray) { Py_DECREF(sendArray); }
    return PyArray_Return(scanArray);
  fail:
    if (is_new_object && sendArray) { Py_DECREF(sendArray); }
    Py_XDECREF(scanArray);
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
%ignore Teuchos::Comm::isend;
%ignore Teuchos::Comm::ireceive;
%ignore Teuchos::Comm::wait;
%ignore Teuchos::Comm::waitAll;
%ignore Teuchos::Comm::readySend;
%ignore Teuchos::broadcast(const Comm&, const int, const ArrayView&);
%ignore Teuchos::wait;
%ignore Teuchos::waitAll;
%include "Teuchos_Comm.hpp"
%template(Comm_long) Teuchos::Comm<long>;
%pythoncode
%{
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
%pythoncode
%{
SerialComm = SerialComm_long
%}

//////////////////////////////////
// Teuchos::Comm Helper support //
//////////////////////////////////
%rename(reductionTypeToString) Teuchos::toString;
%include "Teuchos_CommHelpers.hpp"
%template(rank_long   ) Teuchos::rank<long>;
%template(size_long   ) Teuchos::size<long>;
%template(barrier_long) Teuchos::barrier<long>;
%pythoncode
%{
rank    = rank_long
size    = size_long
barrier = barrier_long

def broadcast(comm, rootRank, buffer):
  """
  broadcast(Comm comm, int rootRank, numpy.ndarray buffer)

  Broadcast the contents of buffer from processor rootRank to all of
  the other processors.  Argument buffer must be a numpy array, so
  that the broadcast can be performed in-place.  Its scalar data type
  can be any numerical type supported by numpy.
  """
  comm.broadcast(rootRank, buffer)

def gatherAll(comm, buffer):
  """
  gatherAll(Comm comm, buffer) -> numpy.ndarray

  Gather the contents of buffer to all of the processors.  Argument
  buffer can be a numpy array or any sequence that can be converted to
  a numpy array.  Its scalar data type can be any numerical type
  supported by numpy.  The return argument is a numpy array of the
  same type.
  """
  return comm.gatherAll(buffer)

def reduceAll(comm, reductOp, buffer):
  """
  reduceAll(Comm comm, EReductionType reductOp, buffer) -> numpy.ndarray

  Reduce the contents of buffer according to the operation designated
  by reductOp on all of the processors.  Argument reductOp can be
  Teuchos.REDUCE_SUM, Teuchos.REDUCE_MAX, or Teuchos.REDUCE_MIN.
  Argument buffer can be a numpy array or any sequence that can be
  converted to a numpy array.  Its scalar data type can be any
  numerical type supported by numpy.  The return argument is a numpy
  array of the same type.
  """
  return comm.reduceAll(reductOp, buffer)

def scan(comm, reductOp, buffer):
  """
  scan(Comm comm, EReductionType reductOp, buffer) -> numpy.ndarray

  Return the scan of the contents of buffer according to the operation
  designated by reductOp on each of the processors.  Argument reductOp
  can be Teuchos.REDUCE_SUM, Teuchos.REDUCE_MAX, or
  Teuchos.REDUCE_MIN.  Argument buffer can be a numpy array or any
  sequence that can be converted to a numpy array.  Its scalar data
  type can be any numerical type supported by numpy.  The return
  argument is a numpy array of the same type.
  """
  return comm.scan(reductOp, buffer)
%}

////////////////////////////////////////////////////////////////////////////////
// The following code is implemented if HAVE_MPI is defined
////////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_MPI
%{
#include "mpi.h"
#include "Teuchos_OpaqueWrapper.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
%}

%inline
{
  PyObject* Teuchos_MPI_Init_Argv(PyObject *args)
  {
    // Check if MPI is already initialized
    int ierr = 0;
    MPI_Initialized(&ierr);
    if (ierr) Py_RETURN_FALSE;

    // Reconstruct the command-line arguments
    int    argc  = 0;
    char **argv  = 0;
    if (!PySequence_Check(args))
    {
      PyErr_SetString(PyExc_TypeError, "Init_Argv argument must be a sequence");
      goto fail;
    }
    argc = PySequence_Size(args);
    argv = new char*[argc+1];
    for (int i=0; i<argc; ++i)
    {
      PyObject * item = PySequence_GetItem(args, i);
      if (!PyString_Check(item))
      {
	PyErr_SetString(PyExc_TypeError, "Init_Argv argument list contains non-string");
	goto fail;
      }
      argv[i] = PyString_AsString(item);
    }
    argv[argc] = NULL; //Lam 7.0 requires last arg to be NULL

    //Initialize MPI
    ierr = MPI_Init(&argc, &argv);
    if (ierr)
    {
      PyErr_Format(PyExc_RuntimeError, "MPI initialization error %d", ierr);
      goto fail;
    }
    delete [] argv;
    Py_RETURN_TRUE;
  fail:
    if (argv) delete [] argv;
    return NULL;
  }

  PyObject* Teuchos_MPI_Finalize()
  {
    // Check if MPI has already been finalized
    int ierr = 0;
    MPI_Finalized(&ierr);
    if (ierr) return Py_BuildValue("");

    // Finalize MPI
    ierr = MPI_Finalize();
    if (ierr)
    {
      PyErr_Format(PyExc_RuntimeError, "MPI finalization error %d", ierr);
      return NULL;
    }
    return Py_BuildValue("");
  }
}

// Add python code to call MPI_Init() if appropriate.  If Init_Argv()
// returns True, then MPI_Init() was not called before and Epetra is
// responsible for calling MPI_Finalize(), via the atexit module.
%pythoncode
%{
# Call MPI_Init if appropriate
import sys
calledMpiInit = Teuchos_MPI_Init_Argv(sys.argv)

# Arrange for MPI_Finalize to be called at exit, if appropriate
if calledMpiInit:
    import atexit
    atexit.register(Teuchos_MPI_Finalize)
%}

//////////////////////////////
// Teuchos::MpiComm support //
//////////////////////////////
%extend Teuchos::MpiComm
{
  MpiComm()
  {
    return new Teuchos::MpiComm<Ordinal>
      (Teuchos::opaqueWrapper((MPI_Comm)MPI_COMM_WORLD));
  }
}
%ignore Teuchos::MpiComm::MpiComm;
%ignore Teuchos::MpiComm::getRawMpiComm;
%ignore Teuchos::MpiComm::broadcast;
%ignore Teuchos::MpiComm::gatherAll;
%ignore Teuchos::MpiComm::reduceAll;
%ignore Teuchos::MpiComm::scan;
%include "Teuchos_DefaultMpiComm.hpp"
%template(MpiComm_long) Teuchos::MpiComm<long>;
%pythoncode
%{
MpiComm = MpiComm_long
%}

///////////////////////////////////////////
// Teuchos.DefaultComm support under MPI //
///////////////////////////////////////////
%pythoncode
%{
class DefaultComm:
    "Encapsulate the default global communicator"
    __defaultComm = MpiComm()
    @classmethod
    def getComm(cls):
        "Return the default global communicator"
        return cls.__defaultComm
%}

////////////////////////////////////////////////////////////////////////////////
// The following code is implemented if HAVE_MPI is not defined
////////////////////////////////////////////////////////////////////////////////

#else

/////////////////////////////////////////////
// Teuchos.DefaultComm support without MPI //
/////////////////////////////////////////////
%pythoncode
%{
class DefaultComm:
    "Encapsulate the default global communicator"
    __defaultComm = SerialComm()
    @classmethod
    def getComm(cls):
        "Return the default global communicator"
        return cls.__defaultComm
%}

#endif
