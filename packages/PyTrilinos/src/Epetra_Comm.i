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
// Epetra include files
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Distributor.h"
#include "Epetra_SerialDistributor.h"

// PyTrilinos include files
#include "PyTrilinos_config.h"
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

// General ignore directive
%ignore *::operator=;

/////////////////////////////////////////////////////////
// Teuchos::RCP<> support for all classes in this file //
/////////////////////////////////////////////////////////
%teuchos_rcp(Epetra_BlockMap)
%teuchos_rcp(Epetra_Comm)
%teuchos_rcp(Epetra_SerialComm)
%teuchos_rcp(Epetra_Distributor)
%teuchos_rcp(Epetra_SerialDistributor)
#ifdef HAVE_MPI
%teuchos_rcp(Epetra_MpiComm)
%teuchos_rcp(Epetra_MpiDistributor)
#endif

/////////////////////////
// Epetra_Comm support //
/////////////////////////
%feature("autodoc",
"Broadcast(self, numpy.ndarray myObj, int root)

Argument myObj must be a numpy array, so that the Broadcast can be
performed in-place.  Its scalar data type must be int, long, double or
string.  In C++, this routine has an integer error return code.  In
python, a non-zero return code is converted to an exception.")
Epetra_Comm::Broadcast;
%feature("docstring")
Epetra_Comm::GatherAll
"Argument myObj can be a numpy array or any sequence that can be
converted to a numpy array.  Its scalar data type must be int, long or
double.  The return argument is a numpy array of the same type.  In
C++, this routine has an integer error return code.  In python, a
non-zero return code is converted to an exception."
%feature("docstring")
Epetra_Comm::SumAll
"Argument myObj can be a numpy array or any sequence that can be
converted to a numpy array.  Its scalar data type must be int, long or
double.  The return argument is a numpy array of the same type.  In
C++, this routine has an integer error return code.  In python, a
non-zero return code is converted to an exception."
%feature("docstring")
Epetra_Comm::MaxAll
"Argument myObj can be a numpy array or any sequence that can be
converted to a numpy array.  Its scalar data type must be int, long or
double.  The return argument is a numpy array of the same type.  In
C++, this routine has an integer error return code.  In python, a
non-zero return code is converted to an exception."
%feature("docstring")
Epetra_Comm::MinAll
"Argument myObj can be a numpy array or any sequence that can be
converted to a numpy array.  Its scalar data type must be int, long or
double.  The return argument is a numpy array of the same type.  In
C++, this routine has an integer error return code.  In python, a
non-zero return code is converted to an exception."
%feature("docstring")
Epetra_Comm::ScanSum
"Argument myObj can be a numpy array or any sequence that can be
converted to a numpy array.  Its scalar data type must be int, long or
double.  The return argument is a numpy array of the same type.  In
C++, this routine has an integer error return code.  In python, a
non-zero return code is converted to an exception."
%rename(Comm) Epetra_Comm;
// Several of the Epetra_Comm methods require the same coding pattern
// for their wrappers and can be collapsed into a macro
%fragment("NumPy_Fragments"); // These following macros depend upon all the NumPy fragments
%define %epetra_comm_reduce_method(methodName)
PyObject* methodName(PyObject* partialObj)
{
  int is_new_object, type, count, result;
  PyObject* globalObj = NULL;
  PyArrayObject* partialArray;
  partialArray= obj_to_array_contiguous_allow_conversion(partialObj, NPY_NOTYPE,
							 &is_new_object);
  if (!partialArray) goto fail;
  type      = array_type(partialArray);
  count     = PyArray_SIZE(partialArray);
  globalObj = PyArray_SimpleNew(array_numdims(partialArray),
				array_dimensions(partialArray), type);
  PyArray_FILLWBYTE((PyArrayObject*)globalObj, 0);
  if (type == NPY_INT)
  {
    int* partialVals = (int*) array_data(partialArray);
    int* globalVals  = (int*) array_data(globalObj);
    result = self->methodName(partialVals,globalVals,count);
  }
  else if (type == NPY_LONG)
  {
    long* partialVals = (long*) array_data(partialArray);
    long* globalVals  = (long*) array_data(globalObj);
    result = self->methodName(partialVals,globalVals,count);
  }
  else if (type == NPY_DOUBLE)
  {
    double* partialVals = (double*) array_data(partialArray);
    double* globalVals  = (double*) array_data(globalObj);
    result = self->methodName(partialVals,globalVals,count);
  }
  else
  {
    PyErr_Format(PyExc_TypeError, "Require int, long or double array, got %s array",
		 typecode_string(type));
    goto fail;
  }
  if (result)
  {
    PyErr_Format(PyExc_RuntimeError, "methodName returned error code %d", result);
    goto fail;
  }
  if (is_new_object && partialArray) { Py_DECREF(partialArray); }
  return PyArray_Return((PyArrayObject*)globalObj);

  fail:
  if (is_new_object && partialArray) { Py_DECREF(partialArray); }
  Py_XDECREF(globalObj);
  return NULL;
}
%enddef

// Many of the communicator methods take C arrays as input or output
// arguments.  Here I allow the python user to use numpy arrays
// instead, and for pure input arrays, any python object that can be
// used to construct a numpy array.  Typemaps are not used because
// these methods are overloaded by array type, and the SWIG
// overloading mechanism cannot disambiguate arrays by type.  I only
// extend the base class (Epetra_Comm), which is where I do type
// checking, and rely on Python polymorphism for the derived
// classes.  Also, I do not return an int, but rather raise an
// exception if the routines return a non-zero error code.  Output
// arrays are moved from the argument list to being returned by the
// method.

%extend Epetra_Comm
{
  PyObject* Broadcast(PyObject* myObj, int root)
  {
    int count, type, result;
    PyArrayObject* myArray  = NULL;
    myArray = obj_to_array_no_conversion(myObj, NPY_NOTYPE);
    if (!myArray || !require_contiguous(myArray)) goto fail;
    count = PyArray_SIZE(myArray);
    type  = array_type(myArray);
    if (type == NPY_INT)
    {
      int* myVals = (int*) array_data(myArray);
      result = self->Broadcast(myVals, count, root);
    }
    else if (type == NPY_LONG)
    {
      long* myVals = (long*) array_data(myArray);
      result = self->Broadcast(myVals, count, root);
    }
    else if (type == NPY_DOUBLE)
    {
      double* myVals = (double*) array_data(myArray);
      result = self->Broadcast(myVals, count, root);
    }
    else if (type == NPY_STRING)
    {
      count = PyArray_NBYTES(myArray);
      char* myVals = (char*) array_data(myArray);
      result = self->Broadcast(myVals, count, root);
    }
    else
    {
      PyErr_Format(PyExc_TypeError,
		   "Require an array of int, long, double or string; got %s array",
		   typecode_string(type));
      goto fail;
    }
    if (result)
    {
      PyErr_Format(PyExc_RuntimeError, "Broadcast returned error code %d", result);
      goto fail;
    }
    return Py_BuildValue("");
  fail:
    return NULL;
  }

  PyObject* GatherAll(PyObject* myObj)
  {
    int is_new_object, type, myCount, allND, result;
    PyObject* allObj = NULL;
    PyArrayObject* myArray;
    myArray = obj_to_array_contiguous_allow_conversion(myObj, NPY_NOTYPE,
						       &is_new_object);
    if (!myArray) goto fail;
    type    = array_type(myArray);
    myCount = PyArray_SIZE(myArray);
    allND   = array_numdims(myArray) + 1;
    { // Scope this to make allDims array temporary
      npy_intp * allDims = new npy_intp[allND];
      allDims[0] = self->NumProc();
      for (int i=1; i<allND; ++i) allDims[i] = array_size(myArray,i-1);
      allObj = PyArray_SimpleNew(allND, allDims, type);
      delete [] allDims;
    }
    if (!allObj) goto fail;
    if (type == NPY_INT)
    {
      int* myVals  = (int*) array_data(myArray);
      int* allVals = (int*) array_data(allObj);
      result = self->GatherAll(myVals,allVals,myCount);
    }
    else if (type == NPY_LONG)
    {
      long* myVals  = (long*) array_data(myArray);
      long* allVals = (long*) array_data(allObj);
      result = self->GatherAll(myVals,allVals,myCount);
    }
    else if (type == NPY_DOUBLE)
    {
      double* myVals  = (double*) array_data(myArray);
      double* allVals = (double*) array_data(allObj);
      result = self->GatherAll(myVals,allVals,myCount);
    }
    else
    {
      PyErr_Format(PyExc_TypeError, "Require int, long or double array, got %s array",
		   typecode_string(type));
      goto fail;
    }
    if (result)
    {
      PyErr_Format(PyExc_RuntimeError, "GatherAll returned error code %d", result);
      goto fail;
    }
    if (is_new_object && myArray) { Py_DECREF(myArray); }
    return PyArray_Return((PyArrayObject*)allObj);

  fail:
    if (is_new_object && myArray) { Py_DECREF(myArray); }
    Py_XDECREF(allObj);
    return NULL;
  }

  %epetra_comm_reduce_method(SumAll )
  %epetra_comm_reduce_method(MaxAll )
  %epetra_comm_reduce_method(MinAll )
  %epetra_comm_reduce_method(ScanSum)
}
// Ignore any future declarations of these function names
%ignore Epetra_Comm::Broadcast const;
%ignore Epetra_Comm::GatherAll const;
%ignore Epetra_Comm::SumAll    const;
%ignore Epetra_Comm::MaxAll    const;
%ignore Epetra_Comm::MinAll    const;
%ignore Epetra_Comm::ScanSum   const;
%include "Epetra_Comm.h"

///////////////////////////////
// Epetra_SerialComm support //
///////////////////////////////
%rename(SerialComm) Epetra_SerialComm;
%include "Epetra_SerialComm.h"

////////////////////////////////
// Epetra_Distributor support //
////////////////////////////////
%rename(Distributor) Epetra_Distributor;
%include "Epetra_Distributor.h"

//////////////////////////////////////
// Epetra_SerialDistributor support //
//////////////////////////////////////
%rename(SerialDistributor) Epetra_Distributor;
%include "Epetra_SerialDistributor.h"

////////////////////////////////////////////////////////////////////////////////
// The following code is implemented if HAVE_MPI is defined
////////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_MPI
// %{
// #include "mpi.h"
// #include "Epetra_MpiComm.h"
// #include "Epetra_MpiDistributor.h"
// %}

%inline
%{
  PyObject* Epetra_MPI_Init_Argv(PyObject *args)
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

  PyObject* Epetra_MPI_Finalize()
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
// Epetra_MPI_Init_Argv() returns True, then MPI_Init() was not called
// before and Epetra is responsible for calling MPI_Finalize(), via
// the atexit module.
%pythoncode
%{

# Call MPI_Init if appropriate
import sys
calledMpiInit = Epetra_MPI_Init_Argv(sys.argv)

# Proceed according to calledMpiInit.  If calledMpiInit is true, then register a
# call to MPI_Finalize() with the atexit module and use the default value for
# mpiCommunicator, equivalent to MPI_COMM_WORLD.  If calledMpiInit is false, try
# to assess what package is responsible for calling MPI_Init(), currently either
# distarray or mpi4py, and extract the appropriate value for mpiCommunicator.
mpiCommunicator = None
if calledMpiInit:
    import atexit
    atexit.register(Epetra_MPI_Finalize)
else:
    if sys.modules.get("distarray.localapi.mpiutils"):
        dlm = sys.modules["distarray.localapi.mpiutils"]
        mpiCommunicator = dlm.get_base_comm()
    elif sys.modules.get("mpi4py.MPI"):
        MPI = sys.modules["mpi4py.MPI"]
        mpiCommunicator = MPI.COMM_WORLD
%}


////////////////////////////
// Epetra_MpiComm support //
////////////////////////////
%extend Epetra_MpiComm
{
#ifdef HAVE_MPI4PY
  Epetra_MpiComm(MPI_Comm mpiComm = MPI_COMM_WORLD)
  {
    Epetra_MpiComm * result = new Epetra_MpiComm(mpiComm);
    return result;
  }
#else
  Epetra_MpiComm(PyObject * dummy = NULL)
  {
    return new Epetra_MpiComm(MPI_COMM_WORLD);
  }
#endif
}
%ignore Epetra_MpiComm::Epetra_MpiComm;
%rename(MpiComm) Epetra_MpiComm;
%include "Epetra_MpiComm.h"

///////////////////////////////////
// Epetra_MpiDistributor support //
///////////////////////////////////
%rename(MpiDistributor) Epetra_MpiDistributor;
%include "Epetra_MpiDistributor.h"

/////////////////////////////////////
// Epetra.PyComm support under MPI //
/////////////////////////////////////
%pythoncode
%{
def PyComm():
    "PyComm() -> Epetra.MpiComm"
    return MpiComm(mpiCommunicator);
%}

#else

////////////////////////////////////////////////////////////////////////////////
// The following code is implemented if HAVE_MPI is not defined
////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////
// Epetra.PyComm support without MPI //
///////////////////////////////////////
%pythoncode
%{
def PyComm():
  "PyComm() -> Epetra.SerialComm"
  return SerialComm();
%}

#endif
