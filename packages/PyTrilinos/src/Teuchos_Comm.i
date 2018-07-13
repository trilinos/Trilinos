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

%{
// System include files
#include <string>

// Teuchos include files
#include "PyTrilinos_Teuchos_Headers.hpp"
using Teuchos::OpaqueWrapper;
%}

// Convey the PyTrilinos configuration to SWIG
%include "PyTrilinos_config.h"

// Handle Mpi4Py, if we have it
#ifdef HAVE_MPI4PY
%{
#include "mpi4py/mpi4py.h"
%}
%include "mpi4py/mpi4py.i"
%mpi4py_typemap(Comm, MPI_Comm);
#endif

// Teuchos Array support
%include "Teuchos_Array.i"

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
%teuchos_rcp(Teuchos::LabeledObject)
%feature("director") Teuchos::LabeledObject;
%include "Teuchos_LabeledObject.hpp"

//////////////////////////////////
// Teuchos::Describable support //
//////////////////////////////////
%teuchos_rcp(Teuchos::Describable)
%feature("director") Teuchos::Describable;
%include "Teuchos_Describable.hpp"

//////////////////////////////////
// Teuchos::ReductionOp support //
//////////////////////////////////
%include "Teuchos_ReductionOp.hpp"

///////////////////////////
// Teuchos::Comm support //
///////////////////////////
%teuchos_rcp(Teuchos::Comm< int >)
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
Teuchos::Comm::scan;
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
%ignore Teuchos::Comm::broadcast(const int rootRank,
                                 const Ordinal bytes,
                                 char buffer[]) const;
%ignore Teuchos::Comm::gatherAll(const Ordinal sendBytes,
                                 const char sendBuffer[],
                                 const Ordinal recvBytes,
                                 char recvBuffer[]) const;
%ignore Teuchos::Comm::reduceAll(const ValueTypeReductionOp<Ordinal,char> &reductOp,
                                 const Ordinal bytes,
                                 const char sendBuffer[],
                                 char globalReducts[]) const;
%ignore Teuchos::Comm::scan(const ValueTypeReductionOp<Ordinal,char> &reductOp,
                            const Ordinal bytes,
                            const char sendBuffer[],
                            char scanReducts[]) const;
%ignore Teuchos::Comm::isend(const ArrayView<const char> &sendBuffer,
                             const int destRank) const;
%ignore Teuchos::Comm::ireceive(const ArrayView<char> &recvBuffer,
                                const int sourceRank) const;
%ignore Teuchos::Comm::wait(const Ptr<RCP<CommRequest<Ordinal> > >& request) const;
%ignore Teuchos::Comm::waitAll(const ArrayView<RCP<CommRequest<Ordinal> > > &requests) const;
%ignore Teuchos::Comm::waitAll(const ArrayView<RCP<CommRequest<Ordinal> > >& requests,
                               const ArrayView<RCP<CommStatus<Ordinal> > >& statuses) const;
%ignore Teuchos::Comm::readySend(const ArrayView<const char> &sendBuffer,
                                 const int destRank) const;
%ignore Teuchos::broadcast(const Comm&,
                           const int,
                           const ArrayView&);
%ignore Teuchos::wait;
%ignore Teuchos::waitAll;
%include "Teuchos_Comm.hpp"
%template(Comm_int) Teuchos::Comm<int>;
%pythoncode
%{
Comm = Comm_int
%}

/////////////////////////////////
// Teuchos::SerialComm support //
/////////////////////////////////
%teuchos_rcp(Teuchos::SerialComm< int >)
%ignore Teuchos::SerialComm::broadcast;
%ignore Teuchos::SerialComm::gatherAll;
%ignore Teuchos::SerialComm::reduceAll;
%ignore Teuchos::SerialComm::scan;
%include "Teuchos_DefaultSerialComm.hpp"
%template(SerialComm_int) Teuchos::SerialComm<int>;
%pythoncode
%{
SerialComm = SerialComm_int
%}

/////////////////////////////////////
// Teuchos::EReductionType support //
/////////////////////////////////////
%rename(reductionTypeToString) Teuchos::toString;
%include "Teuchos_EReductionType.hpp"

//////////////////////////////////
// Teuchos::Comm Helper support //
//////////////////////////////////
%include "Teuchos_CommHelpers.hpp"
%template(rank_int   ) Teuchos::rank<int>;
%template(size_int   ) Teuchos::size<int>;
%template(barrier_int) Teuchos::barrier<int>;
%pythoncode
%{

rank    = rank_int
size    = size_int
barrier = barrier_int

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
#include "Teuchos_DefaultMpiComm.hpp"
%}

%inline
%{
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
      if (!PyUnicode_Check(item) && !PyString_Check(item))
      {
	PyErr_SetString(PyExc_TypeError, "Init_Argv argument list contains non-string");
	goto fail;
      }
#if PY_VERSION_HEX >= 0x03000000
      PyObject * pyBytes = PyUnicode_AsASCIIString(item);
      if (!pyBytes) goto fail;
      char * bytes = PyBytes_AsString(pyBytes);
      argv[i] = new char[strlen(bytes)+1];
      strcpy(bytes,argv[i]);
      Py_DECREF(pyBytes);
#else
      argv[i] = PyString_AsString(item);
#endif
      if (PyErr_Occurred()) goto fail;
    }
    argv[argc] = NULL; //Lam 7.0 requires last arg to be NULL

    //Initialize MPI
    ierr = MPI_Init(&argc, &argv);
    if (ierr)
    {
      PyErr_Format(PyExc_RuntimeError, "MPI initialization error %d", ierr);
      goto fail;
    }
#if PY_VERSION_HEX >= 0x03000000
    for (int i=0; i<argc; ++i) delete [] argv[i];
#endif
    delete [] argv;
    Py_RETURN_TRUE;

  fail:
    if (argv)
    {
#if PY_VERSION_HEX >= 0x03000000
      for (int i=0; i<argc; ++i)
        if (argv[i]) delete [] argv[i];
#endif
      delete [] argv;
    }
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
%}

// Add python code to call MPI_Init() if appropriate.  If
// Teuchos_MPI_Init_Argv() returns True, then MPI_Init() was not
// called before and Epetra is responsible for calling MPI_Finalize(),
// via the atexit module.
%pythoncode
%{

# Call MPI_Init if appropriate
import sys
calledMpiInit = Teuchos_MPI_Init_Argv(sys.argv)

# Proceed according to calledMpiInit.  If calledMpiInit is true, then register a
# call to MPI_Finalize() with the atexit module and use the default value for
# mpiCommunicator, equivalent to MPI_COMM_WORLD.  If calledMpiInit is false, try
# to assess what package is responsible for calling MPI_Init(), currently either
# distarray or mpi4py, and extract the appropriate value for mpiCommunicator.
mpiCommunicator = None
if calledMpiInit:
    import atexit
    atexit.register(Teuchos_MPI_Finalize)
else:
    if sys.modules.get("distarray.localapi.mpiutils"):
        dlm = sys.modules["distarray.localapi.mpiutils"]
        mpiCommunicator = dlm.get_base_comm()
    elif sys.modules.get("mpi4py.MPI"):
        MPI = sys.modules["mpi4py.MPI"]
        mpiCommunicator = MPI.COMM_WORLD
%}

//////////////////////////////
// Teuchos::MpiComm support //
//////////////////////////////
%extend Teuchos::MpiComm
{
#ifdef HAVE_MPI4PY
  MpiComm(MPI_Comm mpiComm = MPI_COMM_WORLD)
  {
    return new Teuchos::MpiComm< Ordinal >
      (Teuchos::opaqueWrapper(mpiComm));
  }
#else
  MpiComm(PyObject * dummy = NULL)
  {
    return new Teuchos::MpiComm< Ordinal >
      (Teuchos::opaqueWrapper(MPI_COMM_WORLD));
  }
#endif
}
%ignore Teuchos::MpiComm::MpiComm(MPI_Comm rawMpiComm);
%ignore Teuchos::MpiComm::MpiComm(const RCP<const OpaqueWrapper<MPI_Comm> >& rawMpiComm);
%ignore Teuchos::MpiComm::MpiComm(const RCP<const OpaqueWrapper<MPI_Comm> >& rawMpiComm,
                                  const int defaultTag);
%ignore Teuchos::MpiComm::MpiComm(const MpiComm<Ordinal>& other);
%ignore Teuchos::MpiComm::getRawMpiComm;
%ignore Teuchos::MpiComm::broadcast;
%ignore Teuchos::MpiComm::gatherAll;
%ignore Teuchos::MpiComm::reduceAll;
%ignore Teuchos::MpiComm::scan;
%teuchos_rcp(Teuchos::MpiComm< int >)
%include "Teuchos_DefaultMpiComm.hpp"
%template(MpiComm_int) Teuchos::MpiComm<int>;
%pythoncode
%{
MpiComm = MpiComm_int
%}

///////////////////////////////////////////
// Teuchos.DefaultComm support under MPI //
///////////////////////////////////////////
%pythoncode
%{

class DefaultComm:
    "Encapsulate the default global communicator"
    __defaultComm = MpiComm(mpiCommunicator)
    @classmethod
    def getComm(cls):
        "Return the default global communicator"
        return cls.__defaultComm
%}

#else

////////////////////////////////////////////////////////////////////////////////
// The following code is implemented if HAVE_MPI is not defined
////////////////////////////////////////////////////////////////////////////////

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
