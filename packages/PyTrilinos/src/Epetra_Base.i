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
// PyTrilinos includes
#include "PythonException.h"
#include "FILEstream.h"

// Epetra includes
#include "Epetra_Version.h"
#include "Epetra_CombineMode.h"
#include "Epetra_DataAccess.h"
#include "Epetra_Object.h"
#include "Epetra_CompObject.h"
#include "Epetra_BLAS.h"
#include "Epetra_LAPACK.h"
#include "Epetra_Flops.h"
#include "Epetra_Time.h"
#include "Epetra_Util.h"
#include "Epetra_MapColoring.h"

// Epetra python exception
char epetraError[13] = "Epetra.Error";
static PyObject * PyExc_EpetraError = PyErr_NewException(epetraError,NULL,NULL);
%}

// NumPy support
%pythoncode
%{
# Much of the Epetra module is compatible with the numpy module
import numpy

# From numpy version 1.0 forward, we want to use the following import syntax for
# user_array.container.  We rename it UserArray for backward compatibility with
# 0.9.x versions of numpy:
try:
    from numpy.lib.user_array import container as UserArray

# If the previous import failed, it is because we are using a numpy version
# prior to 1.0.  So we catch it and try again with different syntax.
except ImportError:

    # There is a bug in UserArray from numpy 0.9.8.  If this is the case, we
    # have our own patched version.
    try:
        from UserArrayFix import UserArray
        
    # If the previous import failed, it is because we are using a version of
    # numpy prior to 0.9.8, such as 0.9.6, which has no bug and therefore has no
    # local patch.  Now we can import using the old syntax:
    except ImportError:
        from numpy.lib.UserArray import UserArray
%}

////////////////////////
// Exception handling //
////////////////////////

// Standard exception handling
%include "exception.i"

// Define the EpetraError exception
%constant PyObject * Error = PyExc_EpetraError;  // This steals the reference

// Director exception handling
%feature("director:except") {
    if ($error != NULL) {
        throw Swig::DirectorMethodException();
    }
}

// General exception handling
%exception
{
  try
  {
    $action
    if (PyErr_Occurred()) SWIG_fail;
  }
  catch(PyTrilinos::PythonException & e)
  {
    e.restore();
    SWIG_fail;
  }
  catch(int errCode)
  {
    PyErr_Format(PyExc_EpetraError, "Error code = %d\nSee stderr for details", errCode);
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch (Swig::DirectorException & e)
  {
    SWIG_fail;
  }
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

// General ignore directives
%ignore *::operator=;                // Not overrideable in python
%ignore *::operator[];               // Replaced with __setitem__ method
%ignore *::operator[] const;         // Replaced with __getitem__ method
%ignore *::UpdateFlops(int) const;   // Use long int version
%ignore *::UpdateFlops(float) const; // Use double version

//////////////////////////////
// Raw data buffer handling //
//////////////////////////////

// Define a macro for converting a method that returns a pointer to
// a 1D array of ints to returning a NumPy array
%define %epetra_intarray1d_output_method(className,methodName,dimMethod)
%ignore className::methodName() const;
%extend className
{
  PyObject * methodName()
  {
    int * result = self->methodName();
    if (result == NULL) return Py_BuildValue("");
    int * data   = NULL;
    npy_intp dims[ ] = { self->dimMethod() };
    PyArray_Descr * dtype = PyArray_DescrFromType(NPY_INT);
    PyObject * returnObj = PyArray_NewFromDescr(&PyArray_Type, dtype, 1, dims, NULL,
						NULL, NPY_FARRAY, NULL);
    if (returnObj == NULL) goto fail;
    data = (int*) array_data(returnObj);
    for (int i=0; i<dims[0]; ++i) data[i] = result[i];
    return PyArray_Return((PyArrayObject*)returnObj);
  fail:
    return NULL;
  }
}
%enddef

// Define a macro for converting a method that returns a pointer to
// a 1D array of doubles to returning a NumPy array
%define %epetra_array1d_output_method(className,methodName,dimMethod)
%ignore className::methodName() const;
%extend className
{
  PyObject * methodName()
  {
    double * result = self->methodName();
    if (result == NULL) return Py_BuildValue("");
    double * data   = NULL;
    npy_intp dims[ ] = { self->dimMethod() };
    PyArray_Descr * dtype = PyArray_DescrFromType(NPY_DOUBLE);
    PyObject * returnObj = PyArray_NewFromDescr(&PyArray_Type, dtype, 1, dims, NULL,
						NULL, NPY_FARRAY, NULL);
    if (returnObj == NULL) goto fail;
    data = (double*) array_data(returnObj);
    for (int i=0; i<dims[0]; ++i) data[i] = result[i];
    return PyArray_Return((PyArrayObject*)returnObj);
  fail:
    return NULL;
  }
}
%enddef

// Define a macro for converting a method that returns a pointer to
// a 2D array of doubles to returning a NumPy array
%define %epetra_array2d_output_method(className,methodName,dimMethod1,dimMethod2)
%ignore className::methodName() const;
%extend className
{
  PyObject * methodName()
  {
    double * result = self->methodName();
    if (result == NULL) return Py_BuildValue("");
    double * data   = NULL;
    npy_intp dims[ ] = { self->dimMethod1(), self->dimMethod2() };
    PyArray_Descr * dtype = PyArray_DescrFromType(NPY_DOUBLE);
    PyObject * returnObj = PyArray_NewFromDescr(&PyArray_Type, dtype, 2, dims, NULL,
						NULL, NPY_FARRAY, NULL);
    if (returnObj == NULL) goto fail;
    data = (double*) array_data(returnObj);
    for (int i=0; i<dims[0]*dims[1]; ++i) data[i] = result[i];
    return PyArray_Return((PyArrayObject*)returnObj);
  fail:
    return NULL;
  }
}
%enddef

// Define macro for typemaps that convert arguments from
// Epetra_*Matrix or Epetra_*Vector to the corresponding
// Epetra_NumPy*Matrix or Epetra_NumPy*Vector.  There is additional
// magic in the python code to convert the Epetra_NumPy*Matrix or
// Epetra_NumPy*Vector to to an Epetra.*Matrix or Epetra.*Vector.
// Since all of these classes potentially represent large data
// buffers, we want efficient memory management and so store them
// internally with Teuchos::RCP<>.
#ifdef HAVE_TEUCHOS
%define %teuchos_rcp_epetra_numpy_overrides(CONST, CLASS...)

// Output a plain pointer
%typemap(out) CONST Epetra_##CLASS *
{
  Teuchos::RCP< CONST PyTrilinos::Epetra_NumPy##CLASS > *smartresult = 0;
  if ($1)
  {
    CONST PyTrilinos::Epetra_NumPy##CLASS * npa = new CONST PyTrilinos::Epetra_NumPy##CLASS(*$1);
    smartresult = new Teuchos::RCP< CONST PyTrilinos::Epetra_NumPy##CLASS >(npa, bool($owner));
  }
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
                                 $descriptor(Teuchos::RCP< PyTrilinos::Epetra_NumPy##CLASS > *),
				 $owner | SWIG_POINTER_OWN));
}

// Output a plain reference
%apply (CONST Epetra_##CLASS *) {CONST Epetra_##CLASS &}

// Input/output of a reference to a pointer
%typemap(in,numinputs=0) Epetra_##CLASS *& (Epetra_##CLASS * _object)
{
  $1 = &_object;
}
%typemap(argout) Epetra_##CLASS *&
{
  PyTrilinos::Epetra_NumPy##CLASS * npa$argnum = new PyTrilinos::Epetra_NumPy##CLASS(**$1);
  Teuchos::RCP< PyTrilinos::Epetra_NumPy##CLASS > *smartobj$argnum =
    new Teuchos::RCP< PyTrilinos::Epetra_NumPy##CLASS >(npa$argnum);
  PyObject * obj = SWIG_NewPointerObj((void*)smartobj$argnum,
			   $descriptor(Teuchos::RCP< PyTrilinos::Epetra_NumPy##CLASS > *),
				      SWIG_POINTER_OWN);
  $result = SWIG_Python_AppendOutput($result,obj);
}

// Director input of a plain reference
%typemap(directorin) CONST Epetra_##CLASS &
%{
  Teuchos::RCP< CONST PyTrilinos::Epetra_NumPy##CLASS > *temp$argnum = new
    Teuchos::RCP< CONST PyTrilinos::Epetra_NumPy##CLASS >
    (new PyTrilinos::Epetra_NumPy##CLASS(View,$1_name), false);
  $input = SWIG_NewPointerObj((void*)temp$argnum,
			      $descriptor(Teuchos::RCP< PyTrilinos::Epetra_NumPy##CLASS > *),
			      SWIG_POINTER_OWN);
%}

%enddef

#else
// HAVE_TEUCHOS is not defined here, so provide %typemaps that assume
// internal storage is via raw pointers.

%define %teuchos_rcp_epetra_numpy_overrides(CONST, CLASS...)

// Output a plain pointer
%typemap(out) CONST Epetra_##CLASS *
{
  if ($1 == NULL) $result = Py_BuildValue("");
  else
  {
    PyTrilinos::Epetra_NumPy##CLASS * npa = new PyTrilinos::Epetra_NumPy##CLASS(*$1);
    $result = SWIG_NewPointerObj(npa, $descriptor(PyTrilinos::Epetra_NumPy##CLASS*), $owner);
  }
}

// Output a plain reference
%apply (CONST Epetra_##CLASS *) {CONST Epetra_##CLASS &}

// Input/output of a reference to a pointer
%typemap(in,numinputs=0) Epetra_##CLASS *& (Epetra_##CLASS * _object)
{
  $1 = &_object;
}
%typemap(argout) Epetra_##CLASS *&
{
  PyTrilinos::Epetra_NumPy##CLASS * npa$argnum = new PyTrilinos::Epetra_NumPy##CLASS(**$1);
  PyObject * obj = SWIG_NewPointerObj((void*)(npa$argnum),
				      $descriptor(PyTrilinos::Epetra_NumPy##CLASS*),
				      SWIG_POINTER_OWN);
  $result = SWIG_Python_AppendOutput($result,obj);
}

// Director input of a plain reference
%typemap(directorin) CONST Epetra_##CLASS &
%{
  PyTrilinos::Epetra_NumPy##CLASS *npa$argnum = new PyTrilinos::Epetra_NumPy##CLASS(View,$1_name);
  $input = SWIG_NewPointerObj(npa$argnum,
			      $descriptor(PyTrilinos::Epetra_NumPy##CLASS*),
			      SWIG_POINTER_OWN);
%}

%enddef
#endif

// Use the %teuchos_rcp_epetra_numpy_overrides macro to define the
// %teuchos_rcp_epetra_numpy macro
#define EMPTYHACK
%define %teuchos_rcp_epetra_numpy(CLASS)
  %teuchos_rcp(Epetra_##CLASS)
  %teuchos_rcp(PyTrilinos::Epetra_NumPy##CLASS)
  %teuchos_rcp_epetra_numpy_overrides(EMPTYHACK, CLASS)
  %teuchos_rcp_epetra_numpy_overrides(const, CLASS)
%enddef

// Define macros for typemaps that convert a reference to a pointer to
// an object, into a python return argument (which might be placed into a
// tuple, if there are more than one).
#ifdef HAVE_TEUCHOS

%define %teuchos_rcp_epetra_argout(ClassName)
%typemap(in,numinputs=0) ClassName *& (ClassName * _object)
{
  $1 = &_object;
}
%typemap(argout) ClassName *&
{
  Teuchos::RCP< ClassName > *smartresult = new Teuchos::RCP< ClassName >(*$1);
  PyObject * obj = SWIG_NewPointerObj(%as_voidptr(smartresult),
				      $descriptor(Teuchos::RCP< ClassName >*), SWIG_POINTER_OWN);
  $result = SWIG_Python_AppendOutput($result,obj);
}
%enddef

#else

%define %teuchos_rcp_epetra_argout(ClassName)
%typemap(in,numinputs=0) ClassName *& (ClassName * _object)
{
  $1 = &_object;
}
%typemap(argout) ClassName *&
{
  PyObject * obj = SWIG_NewPointerObj(%as_voidptr($1), $descriptor(ClassName *), SWIG_POINTER_OWN);
  $result = SWIG_Python_AppendOutput($result,obj);
}
%enddef

#endif

////////////////////////////
// Epetra_Version support //
////////////////////////////
%include "Epetra_Version.h"
%pythoncode
%{
Version = Epetra_Version
__version__ = Version().split()[2]
%}

////////////////////////////////
// Epetra_CombineMode support //
////////////////////////////////
%include "Epetra_CombineMode.h"

///////////////////////////////
// Epetra_DataAccess support //
///////////////////////////////
%include "Epetra_DataAccess.h"

///////////////////////////
// Epetra_Object support //
///////////////////////////
%feature("docstring")
Epetra_Object
"The base Epetra class.
    
The Epetra_Object class provides capabilities common to all Epetra
objects, such as a label that identifies an object instance, constant
definitions, enum types.  In C++, it supports a ``Print()`` method
that takes an output stream as an argument.  In the python
implementation for this and all derived classes, this method takes an
optional file object argument whose default value is standard out."
%feature("docstring")
Epetra_Object::__str__
"Returns the results of ``Print()`` in a string, so that
the python ``print`` command will work on ``Epetra`` objects.  The
``Print()`` methods are designed to run correctly in parallel, so do
not execute ``print`` on an Epetra object conditionally on the
processor number.  For example, do not do

  ``if comm.MyPID() == 0: print epetra_obj``

or it will hang your code."
%teuchos_rcp(Epetra_Object)
%rename(Object) Epetra_Object;
%extend Epetra_Object
{
  // The __str__() method is used by the python str() operator on any
  // object given to the python print command.
  PyObject * __str__()
  {
    std::ostringstream os;
    self->Print(os);               // Put the output in os
    std::string s = os.str();      // Extract the string from os
    Py_ssize_t last = s.length();  // Get the last index
    if (s.substr(last) == "\n")
      last-=1;                     // Ignore any trailing newline
    return PyString_FromStringAndSize(s.c_str(),last);  // Return the string as a PyObject
  }
  // The ClassName::Print(ostream) method will be ignored and replaced
  // by a MyPrint() (renamed Print()) method here that takes a python
  // file object as its optional argument.  If no argument is given,
  // then output is to standard out.
  void Print(PyObject*pf=NULL) const
  {
    if (pf == NULL)
    {
      self->Print(std::cout);
    }
    else
    {
      if (!PyFile_Check(pf))
      {
	PyErr_SetString(PyExc_IOError, "Print() method expects file object");
      }
      else
      {
	std::FILE * f = PyFile_AsFile(pf);
	PyTrilinos::FILEstream buffer(f);
	std::ostream os(&buffer);
	self->Print(os);
	os.flush();
      }
    }
  }
}
// Ignore the Print() method for all classes.  Only the above Print()
// method will be implemented
%ignore *::Print;
%ignore operator<<(ostream &, const Epetra_Object &); // From python, use __str__
%include "Epetra_Object.h"

///////////////////////////////
// Epetra_CompObject support //
///////////////////////////////
%teuchos_rcp(Epetra_CompObject)
%rename(CompObject) Epetra_CompObject;
%include "Epetra_CompObject.h"

/////////////////////////
// Epetra_BLAS support //
/////////////////////////
// I used to %import here, but newer versions of swig raise a bunch of
// warnings for doing this.  Now I use %include, coupled with a bunch
// of %ignores, because I want a simple python base class without the
// C-style BLAS interface.
%teuchos_rcp(Epetra_BLAS)
%rename(BLAS) Epetra_BLAS;
%ignore Epetra_BLAS::ASUM;
%ignore Epetra_BLAS::AXPY;
%ignore Epetra_BLAS::COPY;
%ignore Epetra_BLAS::DOT;
%ignore Epetra_BLAS::GEMM;
%ignore Epetra_BLAS::GEMV;
%ignore Epetra_BLAS::IAMAX;
%ignore Epetra_BLAS::NRM2;
%ignore Epetra_BLAS::SCAL;
%ignore Epetra_BLAS::SYMM;
%ignore Epetra_BLAS::TRMM;
%include "Epetra_BLAS.h"

///////////////////////////
// Epetra_LAPACK support //
///////////////////////////
// I used to %import here, but newer versions of swig raise a bunch of
// warnings for doing this.  Now I use %include, coupled with a bunch
// of %ignores, because I want a simple python base class without the
// C-style LAPACK interface.
%rename(LAPACK) Epetra_LAPACK;
%ignore Epetra_LAPACK::GECON;
%ignore Epetra_LAPACK::GEEQU;
%ignore Epetra_LAPACK::GEEV;
%ignore Epetra_LAPACK::GEEVX;
%ignore Epetra_LAPACK::GEHRD;
%ignore Epetra_LAPACK::GELS;
%ignore Epetra_LAPACK::GEQRF;
%ignore Epetra_LAPACK::GERFS;
%ignore Epetra_LAPACK::GESDD;
%ignore Epetra_LAPACK::GESV;
%ignore Epetra_LAPACK::GESVD;
%ignore Epetra_LAPACK::GESVX;
%ignore Epetra_LAPACK::GETRF;
%ignore Epetra_LAPACK::GETRI;
%ignore Epetra_LAPACK::GETRS;
%ignore Epetra_LAPACK::GGEV;
%ignore Epetra_LAPACK::GGLSE;
%ignore Epetra_LAPACK::GGSVD;
%ignore Epetra_LAPACK::HSEQR;
%ignore Epetra_LAPACK::LAMCH;
%ignore Epetra_LAPACK::LARFT;
%ignore Epetra_LAPACK::ORGHR;
%ignore Epetra_LAPACK::ORGQR;
%ignore Epetra_LAPACK::ORMHR;
%ignore Epetra_LAPACK::POCON;
%ignore Epetra_LAPACK::POEQU;
%ignore Epetra_LAPACK::PORFS;
%ignore Epetra_LAPACK::POSV;
%ignore Epetra_LAPACK::POSVX;
%ignore Epetra_LAPACK::POTRF;
%ignore Epetra_LAPACK::POTRI;
%ignore Epetra_LAPACK::POTRS;
%ignore Epetra_LAPACK::SPEV;
%ignore Epetra_LAPACK::SPGV;
%ignore Epetra_LAPACK::SYEV;
%ignore Epetra_LAPACK::SYEVD;
%ignore Epetra_LAPACK::SYEVR;
%ignore Epetra_LAPACK::SYEVX;
%ignore Epetra_LAPACK::SYGV;
%ignore Epetra_LAPACK::SYGVX;
%ignore Epetra_LAPACK::TREVC;
%ignore Epetra_LAPACK::TREXC;
%include "Epetra_LAPACK.h"

//////////////////////////
// Epetra_Flops support //
//////////////////////////
%rename(FLOPS) Epetra_Flops;
%include "Epetra_Flops.h"

/////////////////////////
// Epetra_Time support //
/////////////////////////
%teuchos_rcp(Epetra_Time)
%rename(Time) Epetra_Time;
%include "Epetra_Time.h"

