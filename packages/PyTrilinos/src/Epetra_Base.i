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
// PyTrilinos include files
#include "PyTrilinos_config.h"
#include "PyTrilinos_PythonException.hpp"
#include "PyTrilinos_FILEstream.hpp"

// Epetra include files
#include "Epetra_Version.h"
#include "PyTrilinos_Epetra_Headers.hpp"

// Epetra python exception
char epetraError[13] = "Epetra.Error";
static PyObject * PyExc_EpetraError = PyErr_NewException(epetraError,NULL,NULL);
%}

// NumPy support
%pythoncode
%{
# Much of the Epetra module is compatible with the numpy module
import numpy
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
  catch (Swig::DirectorException & e)
  {
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
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
						NULL, NPY_ARRAY_FARRAY, NULL);
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
						NULL, NPY_ARRAY_FARRAY, NULL);
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
						NULL, NPY_ARRAY_FARRAY, NULL);
    if (returnObj == NULL) goto fail;
    data = (double*) array_data(returnObj);
    for (int i=0; i<dims[0]*dims[1]; ++i) data[i] = result[i];
    return PyArray_Return((PyArrayObject*)returnObj);
  fail:
    return NULL;
  }
}
%enddef

// // Define macro for typemaps that convert arguments from
// // Epetra_*Matrix or Epetra_*Vector to the corresponding
// // Epetra_NumPy*Matrix or Epetra_NumPy*Vector.  There is additional
// // magic in the python code to convert the Epetra_NumPy*Matrix or
// // Epetra_NumPy*Vector to to an Epetra.*Matrix or Epetra.*Vector.
// // Since all of these classes potentially represent large data
// // buffers, we want efficient memory management and so store them
// // internally with Teuchos::RCP<>.
// %define %teuchos_rcp_epetra_numpy_overrides(CONST, CLASS...)

// // Output a plain pointer
// %typemap(out) CONST Epetra_##CLASS *
// {
//   Teuchos::RCP< CONST PyTrilinos::Epetra_NumPy##CLASS > *smartresult = 0;
//   if ($1)
//   {
//     CONST PyTrilinos::Epetra_NumPy##CLASS * npa = new CONST PyTrilinos::Epetra_NumPy##CLASS(*$1);
//     smartresult = new Teuchos::RCP< CONST PyTrilinos::Epetra_NumPy##CLASS >(npa, bool($owner));
//   }
//   %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
//                                  $descriptor(Teuchos::RCP< PyTrilinos::Epetra_NumPy##CLASS > *),
// 				 $owner | SWIG_POINTER_OWN));
// }

// // Output a plain reference
// %apply (CONST Epetra_##CLASS *) {CONST Epetra_##CLASS &}

// // Input/output of a reference to a pointer
// %typemap(in,numinputs=0) Epetra_##CLASS *& (Epetra_##CLASS * _object)
// {
//   $1 = &_object;
// }
// %typemap(argout) Epetra_##CLASS *&
// {
//   PyTrilinos::Epetra_NumPy##CLASS * npa$argnum = new PyTrilinos::Epetra_NumPy##CLASS(**$1);
//   Teuchos::RCP< PyTrilinos::Epetra_NumPy##CLASS > *smartobj$argnum =
//     new Teuchos::RCP< PyTrilinos::Epetra_NumPy##CLASS >(npa$argnum);
//   PyObject * obj = SWIG_NewPointerObj((void*)smartobj$argnum,
// 			   $descriptor(Teuchos::RCP< PyTrilinos::Epetra_NumPy##CLASS > *),
// 				      SWIG_POINTER_OWN);
//   $result = SWIG_Python_AppendOutput($result,obj);
// }

// // Director input of a plain reference
// %typemap(directorin) CONST Epetra_##CLASS &
// %{
//   Teuchos::RCP< CONST PyTrilinos::Epetra_NumPy##CLASS > *temp$argnum = new
//     Teuchos::RCP< CONST PyTrilinos::Epetra_NumPy##CLASS >
//     (new PyTrilinos::Epetra_NumPy##CLASS(View,$1_name), false);
//   $input = SWIG_NewPointerObj((void*)temp$argnum,
// 			      $descriptor(Teuchos::RCP< PyTrilinos::Epetra_NumPy##CLASS > *),
// 			      SWIG_POINTER_OWN);
// %}

// %enddef

// Use the %teuchos_rcp_epetra_numpy_overrides macro to define the
// %teuchos_rcp_epetra_numpy macro
// %define %teuchos_rcp_epetra_numpy(CLASS)
//   %teuchos_rcp(Epetra_##CLASS)
//   %teuchos_rcp(PyTrilinos::Epetra_NumPy##CLASS)
//   %teuchos_rcp_epetra_numpy_overrides(SWIGEMPTYHACK, CLASS)
//   %teuchos_rcp_epetra_numpy_overrides(const        , CLASS)
// %enddef

// Use the %teuchos_rcp_dap and %teuchos_rcp_typemaps_overrides macros
// to define the %teuchos_rcp_epetra_dap macro
%define %teuchos_rcp_epetra_dap(CLASS...)
  %feature("smartptr", noblock=1) Epetra_##CLASS
  {
    Teuchos::RCP< Epetra_##CLASS >
  }
  %teuchos_rcp_dap(PyTrilinos::convertPythonToEpetra##CLASS, Epetra_##CLASS)
  %teuchos_rcp_typemaps_overrides(SWIGEMPTYHACK, Epetra_##CLASS)
  %teuchos_rcp_typemaps_overrides(const        , Epetra_##CLASS)
%enddef

// Define macros for typemaps that convert a reference to a pointer to
// an object, into a python return argument (which might be placed into a
// tuple, if there are more than one).
%define %teuchos_rcp_epetra_argout(ClassName)
%typemap(in,numinputs=0) ClassName *& (ClassName * _object)
{
  $1 = &_object;
}
%typemap(argout) ClassName *&
{
  Teuchos::RCP< ClassName > *smartresult = new Teuchos::RCP< ClassName >(*$1);
  PyObject * obj = SWIG_NewPointerObj(%as_voidptr(smartresult),
				      $descriptor(Teuchos::RCP< ClassName >*),
                                      SWIG_POINTER_OWN);
  $result = SWIG_Python_AppendOutput($result,obj);
}
%enddef

//////////////////////////////
// Python utility functions //
//////////////////////////////
%pythoncode
%{
  def class_array_inplace_op(self, op_str, other):
    in_op = getattr(self.array, "__i"+op_str+"__")
    in_op(other.array)
    return self

  def class_array_math_op(self, op_str, other):
    # Initialize the result by calling the copy constructor
    result = self.__class__(self)
    # Get the equivalent in-place operator for the result
    in_op = getattr(result.array, "__i"+op_str+"__")
    try:
      in_op(other.array)
    except AttributeError:
      in_op(other)
    return result

  def class_array_rmath_op(self, op_str, other):
    # Initialize the result by calling the copy constructor
    result = self.__class__(self)
    indices = (slice(None),) * len(self.array.shape)
    result.array[indices] = other
    in_op = getattr(result.array, "__i"+op_str+"__")
    in_op(self.array)
    return result

  def class_array_add_math_ops(cls, op_str):
    setattr(cls,
            "__i"+op_str+"__",
            lambda self, other: class_array_inplace_op(self, op_str, other))
    setattr(cls,
            "__"+op_str+"__",
            lambda self, other: class_array_math_op(self, op_str, other))
    setattr(cls,
            "__r"+op_str+"__",
            lambda self, other: class_array_rmath_op(self, op_str, other))

  def class_array_add_math(cls):
    class_array_add_math_ops(cls, "add")
    class_array_add_math_ops(cls, "sub")
    class_array_add_math_ops(cls, "mul")
    class_array_add_math_ops(cls, "add")

  def class_array_comp_op(self, op_str, other):
    comp_op = getattr(self.array, "__"+op_str+"__")
    try:
      return comp_op(other.array)
    except AttributeError:
      return comp_op(other)

  def class_array_add_comp_op(cls, op_str):
    setattr(cls,
            "__"+op_str+"__",
            lambda self, other: class_array_comp_op(self, op_str, other))

  def class_array_add_comp(cls):
    class_array_add_comp_op(cls, "lt")
    class_array_add_comp_op(cls, "le")
    class_array_add_comp_op(cls, "eq")
    class_array_add_comp_op(cls, "ne")
    class_array_add_comp_op(cls, "gt")
    class_array_add_comp_op(cls, "ge")
%}


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
  void Print(PyObject * pf=NULL) const
  {
    if (pf == NULL)
    {
      self->Print(std::cout);
    }
    else
    {
      std::ostringstream s;
      self->Print(s);
      if (PyFile_WriteString(s.str().c_str(), pf))
        throw PyTrilinos::PythonException();
    }
  }
}
// Ignore the Print() method for all classes.  Only the above Print()
// method will be implemented
%ignore *::Print;
%ignore operator<<(std::ostream &, const Epetra_Object &); // From python, use __str__
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
%teuchos_rcp(Epetra_LAPACK)
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

// The Epetra_Time constructor takes an Epetra_Comm as its argument,
// so it needs to know how to convert a PyObject to an Epetra_Comm.
// Since the Epetra_Comm classes get wrapped later in Epetra_Comm.i,
// we set up the Teuchos::RCP typemaps for Epetra_Comm here.
%teuchos_rcp(Epetra_Comm)

/////////////////////////
// Epetra_Time support //
/////////////////////////
%teuchos_rcp(Epetra_Time)
%rename(Time) Epetra_Time;
%include "Epetra_Time.h"
