// -*- C++ -*-

%{
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"

PyObject* Init_Argv(PyObject *args) {
  int ierr = 0;
  MPI_Initialized(&ierr);
  if (ierr)
    return Py_BuildValue("");

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

MPI_Comm CommWorld() {
  return(MPI_COMM_WORLD);
}
#else
PyObject* Init_Argv(PyObject *args) {
  return Py_BuildValue("");
}

PyObject* Finalize() {
  return Py_BuildValue("");
}
#endif
%}

// Cretae python interfaces to MPI initialization and finalization
PyObject* Init_Argv(PyObject *args);
PyObject* Finalize();
#ifdef HAVE_MPI
MPI_Comm CommWorld();
#endif

// Ignore directives
%ignore *::Broadcast(int*   ,int    ,int) const;  // These are replaced by %extend below:
%ignore *::Broadcast(long*  ,int    ,int) const;  //   Broadcast(PyObject*,int)
%ignore *::Broadcast(double*,int    ,int) const;
%ignore *::GatherAll(int*   ,int*   ,int) const;  // These are replaced by %extend below:
%ignore *::GatherAll(long*  ,long*  ,int) const;  //   GatherAll(PyObject*)
%ignore *::GatherAll(double*,double*,int) const;
%ignore *::SumAll(   int*   ,int*   ,int) const;  // These are replaced by %extend below:
%ignore *::SumAll(   long*  ,long*  ,int) const;  //   SumAll(PyObject*)
%ignore *::SumAll(   double*,double*,int) const;
%ignore *::MaxAll(   int*   ,int*   ,int) const;  // These are replaced by %extend below:
%ignore *::MaxAll(   long*  ,long*  ,int) const;  //   MaxAll(PyObject*)
%ignore *::MaxAll(   double*,double*,int) const;
%ignore *::MinAll(   int*   ,int*   ,int) const;  // These are replaced by %extend below:
%ignore *::MinAll(   long*  ,long*  ,int) const;  //   MinAll(PyObject*)
%ignore *::MinAll(   double*,double*,int) const;
%ignore *::ScanSum(  int*   ,int*   ,int) const;  // These are replaced by %extend below:
%ignore *::ScanSum(  long*  ,long*  ,int) const;  //   ScanSum(PyObject*)
%ignore *::ScanSum(  double*,double*,int) const;
%ignore Epetra_SerialComm::operator=(const Epetra_SerialComm &);
#ifdef HAVE_MPI
%ignore Epetra_MpiComm::operator=(const Epetra_MpiComm &);
#endif

// Rename directives
%rename(Comm      ) Epetra_Comm;
%rename(SerialComm) Epetra_SerialComm;
#ifdef HAVE_MPI
%rename(MpiComm   ) Epetra_MpiComm;
#endif

// Include the Numeric typemaps and helper functions
%include "numeric.i"

// Typemaps for interfacing Epetra_Comm objects and python Numeric

// Include directives
%include "Epetra_Comm.h"
%include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
%include "Epetra_MpiComm.h"
#endif

/* Extend directives.  Many of the communicator methods take C arrays
 * as input or output arguments.  These extensions allow the python
 * user to use Numeric arrays instead, and for pure input arrays, any
 * python object that can be used to construct a Numeric array.
 * Typemaps are not used because these methods are overloaded by array
 * type, and the SWIG overloading mechanism cannot disambiguate arrays
 * by type.  I only extend the base class (Epetra_Comm), which is
 * where I do type checking, and rely on python's polymorphism for the
 * derived classes.  Also, I do not return an int, but rather raise an
 * exception if the routines return a non-zero error code.  Output
 * arrays are moved from the argument list to being returned by the
 * method.
 */
%extend Epetra_Comm {

  PyObject* Broadcast(PyObject* myObj, int root) {
    int count, type, result;
    PyArrayObject* myArray;
    myArray = obj_to_array_no_conversion(myObj, PyArray_NOTYPE);
    if (!myArray || !require_contiguous(myArray)) goto fail;
    count = PyArray_SIZE(myArray);
    type  = array_type(myArray);
    if (type == PyArray_INT) {
      int* myVals = (int*)myArray->data;
      result = self->Broadcast(myVals,count,root);
    }
    else if (type == PyArray_LONG) {
      long* myVals = (long*)myArray->data;
      result = self->Broadcast(myVals,count,root);
    }
    else if (type == PyArray_DOUBLE) {
      double* myVals = (double*)myArray->data;
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
    Py_INCREF(Py_None);
    return Py_None;
  fail:
    return NULL;
  }

  PyObject* GatherAll(PyObject* myObj) {
    int is_new_object, type, myCount, allND, result;
    PyObject* allObj;
    PyArrayObject* myArray;
    myArray = obj_to_array_contiguous_allow_conversion(myObj, PyArray_NOTYPE,
						       &is_new_object);
    if (!myArray) goto fail;
    type    = array_type(myArray);
    myCount = PyArray_SIZE(myArray);
    allND   = myArray->nd + 1;
    { // Scope this to make allDims array temporary
      int allDims[allND];
      allDims[0] = self->NumProc();
      for (int i=1; i<allND; ++i) allDims[i] = myArray->dimensions[i-1];
      allObj = PyArray_FromDims(allND, allDims, type);
    }
    if (!allObj) goto fail;
    if (type == PyArray_INT) {
      int* myVals  = (int*)myArray->data;
      int* allVals = (int*)((PyArrayObject*)allObj)->data;
      result = self->GatherAll(myVals,allVals,myCount);
    }
    else if (type == PyArray_LONG) {
      long* myVals  = (long*)myArray->data;
      long* allVals = (long*)((PyArrayObject*)allObj)->data;
      result = self->GatherAll(myVals,allVals,myCount);
    }
    else if (type == PyArray_DOUBLE) {
      double* myVals  = (double*)myArray->data;
      double* allVals = (double*)((PyArrayObject*)allObj)->data;
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
    return allObj;
  fail:
    if (is_new_object && myArray) Py_DECREF(myArray);
    return NULL;
  }

  PyObject* SumAll(PyObject* partialObj) {
    int is_new_object, type, count, result;
    PyObject* globalObj;
    PyArrayObject* partialArray;
    partialArray= obj_to_array_contiguous_allow_conversion(partialObj, PyArray_NOTYPE,
							   &is_new_object);
    if (!partialArray) goto fail;
    type      = array_type(partialArray);
    count     = PyArray_SIZE(partialArray);
    globalObj = PyArray_FromDims(partialArray->nd, partialArray->dimensions, type);
    if (type == PyArray_INT) {
      int* partialVals = (int*)partialArray->data;
      int* globalVals  = (int*)((PyArrayObject*)globalObj)->data;
      result = self->SumAll(partialVals,globalVals,count);
    }
    else if (type == PyArray_LONG) {
      long* partialVals = (long*)partialArray->data;
      long* globalVals  = (long*)((PyArrayObject*)globalObj)->data;
      result = self->SumAll(partialVals,globalVals,count);
    }
    else if (type == PyArray_DOUBLE) {
      double* partialVals = (double*)partialArray->data;
      double* globalVals  = (double*)((PyArrayObject*)globalObj)->data;
      result = self->SumAll(partialVals,globalVals,count);
    }
    else {
      PyErr_Format(PyExc_TypeError, "Require int, long or double array, got %s array",
		   typecode_string(type));
      goto fail;
    }
    if (result) {
      PyErr_Format(PyExc_RuntimeError, "SumAll returned error code %d", result);
      goto fail;
    }
    if (is_new_object && partialArray) Py_DECREF(partialArray);
    return globalObj;
  fail:
    if (is_new_object && partialArray) Py_DECREF(partialArray);
    return NULL;
  }

  PyObject* MaxAll(PyObject* partialObj) {
    int is_new_object, type, count, result;
    PyObject* globalObj;
    PyArrayObject* partialArray;
    partialArray = obj_to_array_contiguous_allow_conversion(partialObj, PyArray_NOTYPE,
							    &is_new_object);
    if (!partialArray) goto fail;
    type      = array_type(partialArray);
    count     = PyArray_SIZE(partialArray);
    globalObj = PyArray_FromDims(partialArray->nd, partialArray->dimensions, type);
    if (type == PyArray_INT) {
      int* partialMaxs = (int*)partialArray->data;
      int* globalMaxs  = (int*)((PyArrayObject*)globalObj)->data;
      result = self->MaxAll(partialMaxs,globalMaxs,count);
    }
    else if (type == PyArray_LONG) {
      long* partialMaxs = (long*)partialArray->data;
      long* globalMaxs  = (long*)((PyArrayObject*)globalObj)->data;
      result = self->MaxAll(partialMaxs,globalMaxs,count);
    }
    else if (type == PyArray_DOUBLE) {
      double* partialMaxs = (double*)partialArray->data;
      double* globalMaxs  = (double*)((PyArrayObject*)globalObj)->data;
      result = self->MaxAll(partialMaxs,globalMaxs,count);
    }
    else {
      PyErr_Format(PyExc_TypeError, "Require int, long or double array, got %s array",
		   typecode_string(type));
      goto fail;
    }
    if (result) {
      PyErr_Format(PyExc_RuntimeError, "MaxAll returned error code %d", result);
      goto fail;
    }
    if (is_new_object && partialArray) Py_DECREF(partialArray);
    return globalObj;
  fail:
    if (is_new_object && partialArray) Py_DECREF(partialArray);
    return NULL;
  }

  PyObject* MinAll(PyObject* partialObj) {
    int is_new_object, type, count, result;
    PyObject* globalObj;
    PyArrayObject* partialArray;
    partialArray = obj_to_array_contiguous_allow_conversion(partialObj, PyArray_NOTYPE,
							    &is_new_object);
    if (!partialArray) goto fail;
    type      = array_type(partialArray);
    count     = PyArray_SIZE(partialArray);
    globalObj = PyArray_FromDims(partialArray->nd, partialArray->dimensions, type);
    if (type == PyArray_INT) {
      int* partialMins = (int*)partialArray->data;
      int* globalMins  = (int*)((PyArrayObject*)globalObj)->data;
      result = self->MinAll(partialMins,globalMins,count);
    }
    else if (type == PyArray_LONG) {
      long* partialMins = (long*)partialArray->data;
      long* globalMins  = (long*)((PyArrayObject*)globalObj)->data;
      result = self->MinAll(partialMins,globalMins,count);
    }
    else if (type == PyArray_DOUBLE) {
      double* partialMins = (double*)partialArray->data;
      double* globalMins  = (double*)((PyArrayObject*)globalObj)->data;
      result = self->MinAll(partialMins,globalMins,count);
    }
    else {
      PyErr_Format(PyExc_TypeError, "Require int, long or double array, got %s array",
		   typecode_string(type));
      goto fail;
    }
    if (result) {
      PyErr_Format(PyExc_RuntimeError, "MinAll returned error code %d", result);
      goto fail;
    }
    if (is_new_object && partialArray) Py_DECREF(partialArray);
    return globalObj;
  fail:
    if (is_new_object && partialArray) Py_DECREF(partialArray);
    return NULL;
  }

  PyObject* ScanSum(PyObject* myObj) {
    int is_new_object, type, count, result;
    PyObject* scanObj;
    PyArrayObject* myArray;
    myArray= obj_to_array_contiguous_allow_conversion(myObj, PyArray_NOTYPE,
						      &is_new_object);
    if (!myArray) goto fail;
    type    = array_type(myArray);
    count   = PyArray_SIZE(myArray);
    scanObj = PyArray_FromDims(myArray->nd, myArray->dimensions, type);
    if (type == PyArray_INT) {
      int* myVals   = (int*)myArray->data;
      int* scanSums = (int*)((PyArrayObject*)scanObj)->data;
      result = self->ScanSum(myVals,scanSums,count);
    }
    else if (type == PyArray_LONG) {
      long* myVals   = (long*)myArray->data;
      long* scanSums = (long*)((PyArrayObject*)scanObj)->data;
      result = self->ScanSum(myVals,scanSums,count);
    }
    else if (type == PyArray_DOUBLE) {
      double* myVals   = (double*)myArray->data;
      double* scanSums = (double*)((PyArrayObject*)scanObj)->data;
      result = self->ScanSum(myVals,scanSums,count);
    }
    else {
      PyErr_Format(PyExc_TypeError, "Require int, long or double array, got %s array",
		   typecode_string(type));
      goto fail;
    }
    if (result) {
      PyErr_Format(PyExc_RuntimeError, "SumAll returned error code %d", result);
      goto fail;
    }
    if (is_new_object && myArray) Py_DECREF(myArray);
    return scanObj;
  fail:
    if (is_new_object && myArray) Py_DECREF(myArray);
    return NULL;
  }
}

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

#ifndef HAVE_MPI
%pythoncode %{
def PyComm():
  "PyComm() -> Epetra.SerialComm"
  return SerialComm();
%}
#else
%pythoncode %{
def PyComm():
  return MpiComm(CommWorld());
%}
#endif
