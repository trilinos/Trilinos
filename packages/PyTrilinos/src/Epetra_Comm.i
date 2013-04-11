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

// General ignore directive
%ignore *::operator=;

// Forward declare typemaps for Teuchos::RCP< Epetra_BlockMap >
%teuchos_rcp(Epetra_BlockMap)

/////////////////////////
// Epetra_Comm support //
/////////////////////////
%teuchos_rcp(Epetra_Comm)
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
// checking, and rely on python's polymorphism for the derived
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
%teuchos_rcp(Epetra_SerialComm)
%rename(SerialComm) Epetra_SerialComm;
%include "Epetra_SerialComm.h"

////////////////////////////////
// Epetra_Distributor support //
////////////////////////////////
%teuchos_rcp(Epetra_Distributor)
%rename(Distributor) Epetra_Distributor;
%include "Epetra_Distributor.h"

//////////////////////////////////////
// Epetra_SerialDistributor support //
//////////////////////////////////////
%teuchos_rcp(Epetra_SerialDistributor)
%rename(SerialDistributor) Epetra_Distributor;
%include "Epetra_SerialDistributor.h"

////////////////////////////////////////////////////////////////////////////////
// The following code is implemented if HAVE_MPI is defined
////////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_MPI
%{
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_MpiDistributor.h"
%}

%inline
{
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
}

// Add python code to call MPI_Init() if appropriate.  If Init_Argv()
// returns True, then MPI_Init() was not called before and Epetra is
// responsible for calling MPI_Finalize(), via the atexit module.
%pythoncode
%{
# Call MPI_Init if appropriate
import sys
calledMpiInit = Epetra_MPI_Init_Argv(sys.argv)

# Arrange for MPI_Finalize to be called at exit, if appropriate
if calledMpiInit:
    import atexit
    atexit.register(Epetra_MPI_Finalize)
%}

// CommWorld used to be a function, but I have changed it to a
// constant in order to suppress a memory leak warning we are getting
// under swig 1.3.28.  (And with luck, I have actually fixed a memory
// leak ;-)
%inline
{
  extern const MPI_Comm CommWorld = (MPI_COMM_WORLD);
}

////////////////////////////
// Epetra_MpiComm support //
////////////////////////////
%teuchos_rcp(Epetra_MpiComm)
%rename(MpiComm) Epetra_MpiComm;
%include "Epetra_MpiComm.h"

///////////////////////////////////
// Epetra_MpiDistributor support //
///////////////////////////////////
%teuchos_rcp(Epetra_MpiDistributor)
%rename(MpiDistributor) Epetra_MpiDistributor;
%include "Epetra_MpiDistributor.h"

/////////////////////////////////////
// Epetra.PyComm support under MPI //
/////////////////////////////////////
%pythoncode
%{
def PyComm():
  "PyComm() -> Epetra.MpiComm(MPI_COMM_WORLD)"
  return MpiComm(cvar.CommWorld);
%}

////////////////////////////////////////////////////////////////////////////////
// The following code is implemented if HAVE_MPI is not defined
////////////////////////////////////////////////////////////////////////////////

#else

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
