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
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Distributor.h"
#include "Epetra_SerialDistributor.h"
%}

// Create python interfaces to MPI initialization and finalization
PyObject* Init_Argv(PyObject *args);
PyObject* Finalize();

// General ignore directive
%ignore *::operator=;

/////////////////////////
// Epetra_Comm support //
/////////////////////////
%rename(Comm) Epetra_Comm;
// Several of the Epetra_Comm methods require the same coding pattern
// for their wrappers and can be collapsed into a macro
%define %epetra_comm_reduce_method(methodName)
  PyObject* methodName(PyObject* partialObj) {
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
    if (type == NPY_INT) {
      int* partialVals = (int*) array_data(partialArray);
      int* globalVals  = (int*) array_data(globalObj);
      result = self->methodName(partialVals,globalVals,count);
    }
    else if (type == NPY_LONG) {
      long* partialVals = (long*) array_data(partialArray);
      long* globalVals  = (long*) array_data(globalObj);
      result = self->methodName(partialVals,globalVals,count);
    }
    else if (type == NPY_DOUBLE) {
      double* partialVals = (double*) array_data(partialArray);
      double* globalVals  = (double*) array_data(globalObj);
      result = self->methodName(partialVals,globalVals,count);
    }
    else {
      PyErr_Format(PyExc_TypeError, "Require int, long or double array, got %s array",
		   typecode_string(type));
      goto fail;
    }
    if (result) {
      PyErr_Format(PyExc_RuntimeError, "methodName returned error code %d", result);
      goto fail;
    }
    if (is_new_object && partialArray) Py_DECREF(partialArray);
    return PyArray_Return((PyArrayObject*)globalObj);
  fail:
    if (is_new_object && partialArray) Py_DECREF(partialArray);
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
// checking, and rely on python's polymorphism for the derived
// classes.  Also, I do not return an int, but rather raise an
// exception if the routines return a non-zero error code.  Output
// arrays are moved from the argument list to being returned by the
// method.

%extend Epetra_Comm {

  PyObject* Broadcast(PyObject* myObj, int root) {
    int count, type, result;
    PyArrayObject* myArray = NULL;
    myArray = obj_to_array_no_conversion(myObj, NPY_NOTYPE);
    if (!myArray || !require_contiguous(myArray)) goto fail;
    count = PyArray_SIZE(myArray);
    type  = array_type(myArray);
    if (type == NPY_INT) {
      int* myVals = (int*) array_data(myArray);
      result = self->Broadcast(myVals,count,root);
    }
    else if (type == NPY_LONG) {
      long* myVals = (long*) array_data(myArray);
      result = self->Broadcast(myVals,count,root);
    }
    else if (type == NPY_DOUBLE) {
      double* myVals = (double*) array_data(myArray);
      result = self->Broadcast(myVals,count,root);
    }
    else {
      PyErr_Format(PyExc_TypeError, "Require int, long or double array, got %s array",
		   typecode_string(type));
      goto fail;
    }
    if (result) {
      PyErr_Format(PyExc_RuntimeError, "Broadcast returned error code %d", result);
      goto fail;
    }
    return Py_BuildValue("");
  fail:
    return NULL;
  }

  PyObject* GatherAll(PyObject* myObj) {
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
      intp allDims[allND];
      allDims[0] = self->NumProc();
      for (int i=1; i<allND; ++i) allDims[i] = array_size(myArray,i-1);
      allObj = PyArray_SimpleNew(allND, allDims, type);
    }
    if (!allObj) goto fail;
    if (type == NPY_INT) {
      int* myVals  = (int*) array_data(myArray);
      int* allVals = (int*) array_data(allObj);
      result = self->GatherAll(myVals,allVals,myCount);
    }
    else if (type == NPY_LONG) {
      long* myVals  = (long*) array_data(myArray);
      long* allVals = (long*) array_data(allObj);
      result = self->GatherAll(myVals,allVals,myCount);
    }
    else if (type == NPY_DOUBLE) {
      double* myVals  = (double*) array_data(myArray);
      double* allVals = (double*) array_data(allObj);
      result = self->GatherAll(myVals,allVals,myCount);
    }
    else {
      PyErr_Format(PyExc_TypeError, "Require int, long or double array, got %s array",
		   typecode_string(type));
      goto fail;
    }
    if (result) {
      PyErr_Format(PyExc_RuntimeError, "GatherAll returned error code %d", result);
      goto fail;
    }
    if (is_new_object && myArray) Py_DECREF(myArray);
    return PyArray_Return((PyArrayObject*)allObj);
  fail:
    if (is_new_object && myArray) Py_DECREF(myArray);
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


// Python code.  This will be inserted directly into the python module
%pythoncode %{
# Call MPI_Init if appropriate
import sys
Init_Argv(sys.argv)
del sys

# Arrange for MPI_Finalize to be called at exit, if appropriate
import atexit
atexit.register(Finalize)
%}

////////////////////////////////////////////////////////////////////////////////
// The following code is implemented if HAVE_MPI is defined
////////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_MPI
%{
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_MpiDistributor.h"

PyObject* Init_Argv(PyObject *args) {
  int ierr = 0;
  MPI_Initialized(&ierr);
  if (ierr) return Py_BuildValue("");

  int i, error, myid, size;
  int argc = 0;
  char **argv;
  /* Reconstruct C-commandline */
  argc = PyList_Size(args); //Number of commandline arguments
  argv = (char**) malloc((argc+1)*sizeof(char*));
  for (i=0; i<argc; i++) argv[i] = PyString_AsString(PyList_GetItem(args, i));
  argv[argc] = NULL; //Lam 7.0 requires last arg to be NULL
  error = MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (error) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    PyErr_SetString(PyExc_RuntimeError, "MPI initialization error");
    return NULL;
  }
  return Py_BuildValue("");
}

PyObject* Finalize() {
  int error, myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  // FIXME: add if finalized!
  error = MPI_Finalize();
  if (error) {
    PyErr_SetString(PyExc_RuntimeError, "MPI Finalize error");
    return NULL;
  }
  return Py_BuildValue("");
}
%}

// CommWorld used to be a function, but I have changed it to a
// constant in order to suppress a memory leak warning we are getting
// under swig 1.3.28.  (And with luck, I have actually fixed a memory
// leak ;-)
%inline {extern const MPI_Comm CommWorld = (MPI_COMM_WORLD);}

////////////////////////////
// Epetra_MpiComm support //
////////////////////////////
%rename(MpiComm) Epetra_MpiComm;
%include "Epetra_MpiComm.h"

///////////////////////////////////
// Epetra_MpiDistributor support //
///////////////////////////////////
%rename(MpiDistributor) Epetra_Distributor;
%include "Epetra_MpiDistributor.h"

/////////////////////////////////////
// Epetra.PyComm support under MPI //
/////////////////////////////////////
%pythoncode %{
def PyComm():
  "PyComm() -> Epetra.MpiComm(MPI_COMM_WORLD)"
  return MpiComm(cvar.CommWorld);
%}

////////////////////////////////////////////////////////////////////////////////
// The following code is implemented if HAVE_MPI is not defined
////////////////////////////////////////////////////////////////////////////////

#else

%{
PyObject* Init_Argv(PyObject *args) {
  return Py_BuildValue("");
}

PyObject* Finalize() {
  return Py_BuildValue("");
}
%}

///////////////////////////////////////
// Epetra.PyComm support without MPI //
///////////////////////////////////////
%pythoncode %{
def PyComm():
  "PyComm() -> Epetra.SerialComm"
  return SerialComm();
%}
#endif
