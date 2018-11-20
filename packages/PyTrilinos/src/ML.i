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

%define %ml_docstring
"
PyTrilinos.ML is the python interface to the Trilinos multi-level
preconditioner package ML/MLAPI:

    http://trilinos.sandia.gov/packages/ml

The purpose of ML is to provide multilevel preconditioners to
Trilinos.

ML provides the following user-level classes:

    * BaseObject                - Base class for MLAPI objects
    * CompObject                - FLOP counting base class
    * TimeObject                - Timing base class
    * MultiLevelPreconditioner  - Black-box multilevel preconditioner
    * EpetraBaseOperator        - Base class for interface to Epetra
    * BaseOperator              - Base class for all MLAPI operators
    * Space                     - Defines number of elements and their distribution
    * MultiVector               - MLAPI multivector class
    * Operator                  - MLAPI operator class
    * InverseOperator           - MLAPI inverse operator class
    * PyMatrix                  - Python interface to MLAPI operators

For examples of usage, please consult the following scripts in the
example subdirectory of the PyTrilinos package:

    * exMLAPI.py
    * exMLAPI_Simple.py
    * exMLAPI_AztecOO.py
    * exMLAPI_Iterate.py
    * exMLAPI_PyMatrix.py
    * exMLAPI_Smoother.py
"
%enddef

%module(package   = "PyTrilinos",
	directors = "1",
	autodoc   = "1",
	docstring = %ml_docstring) ML

%{
//  PyTrilinos include files
#include "PyTrilinos_PythonException.hpp"
#include "PyTrilinos_FILEstream.hpp"

// System include files
#include <iostream>
#include <sstream>
#include <vector>

// Configuration include files
#include "PyTrilinos_config.h"

// Teuchos include files
#include "PyTrilinos_Teuchos_Headers.hpp"

// Epetra include files
#ifdef HAVE_PYTRILINOS_EPETRA
#include "PyTrilinos_Epetra_Headers.hpp"

// NumPy include
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
#endif

// IFPACK include files
#ifdef HAVE_IFPACK
#include "PyTrilinos_IFPACK_Headers.hpp"
#endif

// ML include files
#include "PyTrilinos_ML_Headers.hpp"
%}

// Include PyTrilinos configuration
%include "PyTrilinos_config.h"

// Standard exception handling
%include "exception.i"

// Auto-documentation feature
%feature("autodoc", "1");

// Include ML documentation
%include "ML_dox.i"

// SWIG NumPy interface file
%include "numpy.i"

// General ignore directives
%ignore *::operator=;
%ignore *::operator[];

// External Trilinos package imports
%import "Teuchos.i"

#ifdef HAVE_PYTRILINOS_EPETRA
%include "Epetra_RowMatrix_Utils.i"
%ignore Epetra_Version;
%import  "Epetra.i"
#if PY_VERSION_HEX >= 0x030000
%pythoncode
%{
Epetra = PyTrilinos.Epetra
%}
#endif

#ifdef HAVE_IFPACK
%ignore Ifpack_Version;
%import  "IFPACK.i"
#endif
#endif

// Teuchos::RCP handling
%teuchos_rcp(MLAPI::DoubleVector)
%teuchos_rcp(MLAPI::ML_Operator_Box)
#ifdef HAVE_PYTRILINOS_EPETRA
%teuchos_rcp(ML_Epetra::MultiLevelPreconditioner)
#endif

// General exception handling
%feature("director:except")
{
  if ($error != NULL)
  {
    throw Swig::DirectorMethodException();
  }
}

%exception
{
  try
  {
    $action
  }
  catch(PyTrilinos::PythonException &pe)
  {
    pe.restore();
    SWIG_fail;
  }
  catch(Swig::DirectorException &e)
  {
    SWIG_fail;
  }
  catch(int ierr)
  {
    SWIG_exception(SWIG_RuntimeError, "ML exception thrown");
  }
  SWIG_CATCH_STDEXCEPT
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unkown C++ exception");
  }
}

///////////////////////
// ml_config support //
///////////////////////
#undef  PACKAGE_BUGREPORT
%ignore PACKAGE_BUGREPORT;
#undef  PACKAGE_NAME
%ignore PACKAGE_NAME;
#undef  PACKAGE_STRING
%ignore PACKAGE_STRING;
#undef  PACKAGE_TARNAME
%ignore PACKAGE_TARNAME;
#undef  PACKAGE_VERSION
%ignore PACKAGE_VERSION;
%include "ml_config.h"
%include "ml_common.h"

/////////////////////////////
// MLAPI_Workspace support //
/////////////////////////////
%include "MLAPI_Workspace.h"

//////////////////////////////
// MLAPI_BaseObject support //
//////////////////////////////
%rename(Print) MyPrint;
namespace MLAPI
{
%feature("autodoc",
"Print(self, file=None)

Print a description of the object to the given file object.  If no
file object is provided, output is to stdout.")
BaseObject::MyPrint;
%feature("autodoc",
"__str__(self) -> str

Return the description of the object as a string.")
BaseObject::__str__;
%extend BaseObject
{
  PyObject * MyPrint(PyObject * ostream=NULL, const bool verbose=true) const
  {
    if (ostream == NULL)
    {
      self->Print(std::cout, verbose);
    }
    else
    {
      std::ostringstream s;
      self->Print(s, verbose);
      if (PyFile_WriteString(s.str().c_str(), ostream))
        throw PyTrilinos::PythonException();
    }
    return Py_BuildValue("");
  }
  // Define the __str__() method, used by the python str() operator
  // on any object given to the python print command.
  std::string __str__() const
  {
    std::stringstream os;
    self->Print(os);                  // Put the output in os
    std::string s = os.str();         // Extract the string from os
    return s.substr(0,s.length()-1);  // Return the string minus trailing \n
  }
}
}
%ignore *::Print(std::ostream&);
%ignore *::Print(std::ostream&, const bool);
%include "MLAPI_BaseObject.h"

/////////////////////////
// MLAPI_Space support //
/////////////////////////
%apply (int DIM1, int* IN_ARRAY1)
      {(const int NumMyElements, const int* MyGlobalElements)};
%ignore MLAPI::Space::GetRCPMyGlobalElements;
namespace MLAPI
{
%feature("autodoc",
"__init__(self, int numGlobalElements, sequence myGlobalElements) -> Space")
Space::Space(const int NumGlobalElements,
	     const int NumMyElements,
	     const int *MyGlobalElements);
%feature("autodoc",
"GetMyGlobalElements(self) -> Epetra.SerialDenseVector

Return an Epetra.SerialDenseVector of global indexes representing the
elements on this processor.  If the Space is distributed linearly,
this method returns None.")
Space::GetMyGlobalElements;
%extend Space
{
  Epetra_IntSerialDenseVector *  GetMyGlobalElements() const
  {
    Epetra_IntSerialDenseVector * result = self->GetRCPMyGlobalElements().get();
    return result;
  }
}
}
%include "MLAPI_Space.h"
%pythoncode
{
def Space_GetMyGlobalElements(*args):
    """
    GetMyGlobalElements(self) -> Epetra.SerialDenseVector

    Return an Epetra.SerialDenseVector of global indexes representing the
    elements on this processor.
    """
    result = _ML.Space_GetMyGlobalElements(*args)
    if result is None:
        offset = args[0].GetOffset()
        data   = range(offset, offset+args[0].GetNumMyElements())
        result = Epetra.IntSerialDenseVector(data)
    return result
Space.GetMyGlobalElements = Space_GetMyGlobalElements
}

//////////////////////////////
// MLAPI_CompObject support //
//////////////////////////////
%include "MLAPI_CompObject.h"

//////////////////////////////
// MLAPI_TimeObject support //
//////////////////////////////
%include "MLAPI_TimeObject.h"

////////////////////////////////
// MLAPI_BaseOperator support //
////////////////////////////////
// %warnfilter(473)     MLAPI::BaseOperator;
%feature("director") MLAPI::BaseOperator;
%include "MLAPI_BaseOperator.h"

///////////////////////////////
// MLAPI_MultiVector support //
///////////////////////////////
%warnfilter(503) operator<<;  // I can't get %ignore to stop this warning...
%typemap(out) double* GetValues
{
  npy_intp dims = (arg1)->GetMyLength();
  $result = PyArray_SimpleNewFromData(1, &dims, NPY_DOUBLE, (void*) $1);
}
namespace MLAPI
{
%ignore MultiVector::MultiVector(const Space &, double **, const int);
%ignore MultiVector::MultiVector(const Space &, Teuchos::RefCountPtr<DoubleVector>);
%ignore MultiVector::MultiVector(const Space &, std::vector<Teuchos::RefCountPtr<DoubleVector> >);
%ignore MultiVector::operator==;   // Always returns false...
%ignore MultiVector::operator();   // Default behavior can bus error...
%ignore MultiVector::GetValues(const int);
%ignore MultiVector::GetRCPValues;
%ignore MultiVector::SetRCPValues;

%feature("autodoc",
"GetValues(self, int v) -> numpy.ndarray

Returns a 1D numpy array representing the v-th vector in the
MultiVector.")
MultiVector::GetValues;
%feature("autodoc",
"__call__(self, int i, int v=0) -> float

Returns the i-th element of the v-th vector.")
MultiVector::__call__;
%feature("autodoc",
"__add__(self, MultiVector rhs) -> MultiVector

Element-wise addition operator.")
MultiVector::__add__;
%feature("autodoc",
"__sub__(self, MultiVector rhs) -> MultiVector

Element-wise subtraction operator.")
MultiVector::__sub__;
%feature("autodoc",
"__mul__(self, MultiVector rhs) -> float

  Dot-product multiplication operator.
")
MultiVector::__mul__(const MultiVector &) const;
%feature("autodoc",
"__mul__(self, int rhs) -> MultiVector

  Scalar multiplication operators.")
MultiVector::__mul__(const int) const;
%feature("autodoc",
"__rmul__(self, int lhs) -> MultiVector

Scalar 'reverse' multiplication operators.")
MultiVector::__rmul__(const int) const;
%feature("autodoc",
"__setitem__(self, PyObject index, int v, PyObject value)

Argument index can be an integer or slice index into a vector.
Argument v must be an integer specifying the desired vector within the
MultiVector.  The value argument must have a shape that matches the
shape of the index.")
MultiVector::__setitem__;
%feature("autodoc",
"__getitem__(self, PyObject index, int v) -> PyObject

Argument index can be an integer or slice index into a vector.
Argument v must be an integer specifying the desired vector within the
MultiVector.  The returned PyObject will be either a scalar or an
array with the same shape as the index.")
MultiVector::__getitem__;
%feature("autodoc",
"Delete(self, int v)

Deletes the v-th vector within the MultiVector.")
MultiVector::Delete;

%extend MultiVector
{
  MultiVector(const Space & vectorSpace, PyObject * pyValues)
  {
    int is_new_object = 0;
    // Convert pyValues argument to a numpy array
    PyObject * array = (PyObject*)
      obj_to_array_contiguous_allow_conversion(pyValues, NPY_DOUBLE,
					       &is_new_object);
    // Error checking
    if (!array) throw PyTrilinos::PythonException();
    int numDims = array_numdims(array);
    int myVectorLen = PyArray_Size(array);
    if (numDims > 1)
      myVectorLen /= array_size(array, 0);
    if (myVectorLen != vectorSpace.GetNumMyElements())
    {
      PyErr_Format(PyExc_TypeError, "Vector space/numpy array length mismatch:"
		   " %d != %d", vectorSpace.GetNumMyElements(),
		   myVectorLen);
      if (is_new_object) Py_DECREF(array);
      throw PyTrilinos::PythonException();
    }
    // Reshape the array for the case where pyValues represents a 1D
    // array
    if (numDims == 1)
    {
      PyObject * shape = PyTuple_New(2);
      PyTuple_SET_ITEM(shape, 0, PyInt_FromLong(1));
      PyTuple_SET_ITEM(shape, 1, PyInt_FromLong(myVectorLen));
      array = (PyObject*) PyArray_Reshape((PyArrayObject*)array, shape);
      Py_DECREF(shape);
      if (!array) throw PyTrilinos::PythonException();
    }
    // Use the numpy array to construct a MultiVector
    int numVectors = array_size(array, 0);
    MultiVector * result = new MultiVector(vectorSpace, numVectors);
    // Copy the array data into the MultiVector
    double * values = (double*) array_data(array);
    for (int v = 0; v < numVectors; ++v)
    {
      for (int i = 0; i < myVectorLen; ++i)
      {
	result->GetValues(v)[i] = values[v*myVectorLen+i];
      }
    }
    // Release our reference to the array
    if (is_new_object) Py_DECREF(array);
    return result;
  }
  MultiVector __add__(const MultiVector& rhs) const
  {
    if (self->GetVectorSpace() != rhs.GetVectorSpace())
    {
      PyErr_SetString(PyExc_ValueError, "Mismatched Vector Spaces");
      throw PyTrilinos::PythonException();
    }
    MultiVector res = MultiVector(*self);
    for (int v = 0; v < self->GetNumVectors(); ++v)
    {
      double       * resVals = res.GetValues(v);
      const double * rhsVals = rhs.GetValues(v);
      for (int i = 0; i < self->GetMyLength(); ++i)
      {
	resVals[i] += rhsVals[i];
      }
    }
    return res;
  }
  MultiVector __sub__(const MultiVector& rhs) const
  {
    if (self->GetVectorSpace() != rhs.GetVectorSpace())
    {
      PyErr_SetString(PyExc_ValueError, "Mismatched Vector Spaces");
      throw PyTrilinos::PythonException();
    }
    MultiVector res = MultiVector(*self);
    for (int v = 0; v < self->GetNumVectors(); ++v)
    {
      double       * resVals = res.GetValues(v);
      const double * rhsVals = rhs.GetValues(v);
      for (int i = 0; i < self->GetMyLength(); ++i)
      {
	resVals[i] -= rhsVals[i];
      }
    }
    return res;
  }
  double __mul__(const MultiVector& rhs) const
  {
    return self->DotProduct(rhs);
  }
  MultiVector __mul__(const double rhs) const
  {
    MultiVector res = MultiVector(*self);
    for (int v = 0; v < self->GetNumVectors(); ++v)
    {
      double * resVals = res.GetValues(v);
      for (int i = 0; i < self->GetMyLength(); ++i)
      {
	resVals[i] *= rhs;
      }
    }
    return res;
  }
  MultiVector __mul__(const int rhs) const
  {
    return MLAPI_MultiVector___mul____SWIG_1(self, double(rhs));
  }
  MultiVector __rmul__(const double lhs) const
  {
    return MLAPI_MultiVector___mul____SWIG_1(self, lhs);
  }
  MultiVector __rmul__(const int lhs) const
  {
    return MLAPI_MultiVector___mul____SWIG_1(self, double(lhs));
  }
  PyObject * __setitem__(PyObject* args, PyObject * value)
  {
    // Convert args to an index object (could be integer or slice) and
    // an integer index representing the vector number
    PyObject * index  = NULL;
    PyObject * vecNum = NULL;
    int        v      = 0;
    if (!(PyArg_ParseTuple(args, "OO:MLAPI_MultiVector___setitem__",
			   &index, &vecNum) &&
	  SWIG_IsOK(SWIG_AsVal_int(vecNum, &v))))
    {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      throw PyTrilinos::PythonException();
    }
    // Convert the vector specified by v into a numpy array
    npy_intp dim = self->GetMyLength();
    PyObject * array = PyArray_SimpleNewFromData(1,&dim,NPY_DOUBLE,
						 (void*)self->GetValues(v));
    if (!array) throw PyTrilinos::PythonException();
    // Call the __setitem__ method
    char methodName[12] = "__setitem__";
    char format[3]      = "OO";
    PyObject * result = PyObject_CallMethod(array, methodName, format, index, value);
    if (!result) throw PyTrilinos::PythonException();
    // Cleanup and return
    Py_DECREF(array);
    return result;
  }
  PyObject * __getitem__(PyObject* args) const
  {
    // Convert args to an index object (could be integer or slice) and
    // an integer index representing the vector number
    PyObject * index  = NULL;
    PyObject * vecNum = NULL;
    int        v      = 0;
    if (!(PyArg_ParseTuple(args, "OO:MLAPI_MultiVector___setitem__",
			   &index, &vecNum) &&
	  SWIG_IsOK(SWIG_AsVal_int(vecNum, &v))))
    {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      throw PyTrilinos::PythonException();
    }
    // Convert the vector specified by v into a numpy array
    npy_intp dim = self->GetMyLength();
    PyObject * array = PyArray_SimpleNewFromData(1,&dim,NPY_DOUBLE,
						 (void*)self->GetValues(v));
    if (!array) throw PyTrilinos::PythonException();
    // Call the __getitem__ method
    char methodName[12] = "__getitem__";
    char format[2]      = "O";
    PyObject * result = PyObject_CallMethod(array, methodName, format, index);
    if (!result) throw PyTrilinos::PythonException();
    // Cleanup and return
    Py_DECREF(array);
    return result;
  }
  double __call__(const int i, const int v=0) const
  {
    if ((v < 0) or (v > self->GetNumVectors()-1))
    {
      PyErr_Format(PyExc_IndexError, "Invalid vector index %d", v);
      throw PyTrilinos::PythonException();
    }
    if ((i < 0) or (i > self->GetMyLength()-1))
    {
      PyErr_Format(PyExc_IndexError, "Invalid element index %d", i);
      throw PyTrilinos::PythonException();
    }
    return self->operator()(i,v);
  }
}
}
%include "MLAPI_MultiVector.h"

/////////////////////////////////////
// MLAPI_MultiVector_Utils support //
/////////////////////////////////////
%include "MLAPI_MultiVector_Utils.h"

////////////////////////////
// MLAPI_Operator support //
////////////////////////////
%include "MLAPI_Operator.h"
#ifdef HAVE_PYTRILINOS_EPETRA
namespace MLAPI
{
%extend Operator
{
  PyObject* __getitem__(PyObject* args) const
  {
    PyObject * rowObj = NULL;
    PyObject * colObj = NULL;
    int Row           = 0;
    int Col           = 0;
    Epetra_RowMatrix* Matrix = self->GetRCPRowMatrix().get();
    if (SWIG_IsOK(SWIG_AsVal_int(args, &Row)))
    {
      return Epetra_RowMatrix_GetEntries(*Matrix, Row);
    }
    else if (PyArg_ParseTuple(args, "OO:MLAPI_Operator___getitem__",
			      &rowObj, &colObj)      &&
	     SWIG_IsOK(SWIG_AsVal_int(rowObj, &Row)) &&
	     SWIG_IsOK(SWIG_AsVal_int(colObj, &Col))    )
    {
      return Epetra_RowMatrix_GetEntry(*Matrix, Row, Col);
    }
    else
    {
      PyErr_SetString(PyExc_IndexError, "Input argument not supported");
      return NULL;
    }
  }
  MultiVector __mul__(const MultiVector& rhs) const
  {
    return(*self * rhs);
  }
  Operator __add__(const Operator& rhs) const
  {
    return(*self + rhs);
  }
  Operator __sub__(const Operator& rhs) const
  {
    return(*self - rhs);
  }
  Operator __mul__(const Operator& rhs) const
  {
    if (self == &rhs)
      return(*self * Duplicate(rhs));
    else
      return(*self * rhs);
  }
  Operator __mul__(const double rhs) const
  {
    return(*self * rhs);
  }
  Operator __div__(const double rhs) const
  {
    return(*self * (1.0 / rhs));
  }
}
}
#endif

////////////////////////////////
// PyTrilinos_ML_Util support //
////////////////////////////////
#ifdef HAVE_PYTRILINOS_EPETRA
%include "PyTrilinos_ML_Util.hpp"
%extend PyTrilinos::PyMatrix
{
  PyObject * __setitem__(PyObject* args, double val)
  {
    PyObject * rowObj = NULL;
    PyObject * colObj = NULL;
    int Row           = 0;
    int Col           = 0;
    if (!(PyArg_ParseTuple(args, "OO:PyMatrix___setitem__",
			   &rowObj, &colObj)      &&
	  SWIG_IsOK(SWIG_AsVal_int(rowObj, &Row)) &&
	  SWIG_IsOK(SWIG_AsVal_int(colObj, &Col))   ))
    {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return NULL;
    }
    self->SetElement(Row, Col, val);
    return Py_BuildValue("");
  }
  PyObject* __getitem__(PyObject* args)
  {
    PyObject * rowObj = NULL;
    PyObject * colObj = NULL;
    int Row           = 0;
    int Col           = 0;
    if (SWIG_IsOK(SWIG_AsVal_int(args, &Row)))
    {
      return Epetra_RowMatrix_GetEntries(*(self->GetMatrix()), Row);
    }
    else if (PyArg_ParseTuple(args, "OO:PyMatrix___getitem__",
			      &rowObj, &colObj)      &&
	     SWIG_IsOK(SWIG_AsVal_int(rowObj, &Row)) &&
	     SWIG_IsOK(SWIG_AsVal_int(colObj, &Col))    )
    {
      return Epetra_RowMatrix_GetEntry(*(self->GetMatrix()), Row, Col);
    }
    else
    {
      PyErr_SetString(PyExc_IndexError, "Input argument not supported");
      return NULL;
    }
  }
}
#endif

/////////////////////////////////////////
// ml_MultiLevelPreconditioner support //
/////////////////////////////////////////
%include "ml_MultiLevelPreconditioner.h"
#ifdef HAVE_PYTRILINOS_EPETRA
namespace ML_Epetra
{
%extend MultiLevelPreconditioner
{
  // This is for MatrixPortal
  int SetParameterListAndNullSpace(PyObject* obj,
				   const Epetra_MultiVector& NullSpace)
  {
    Teuchos::ParameterList * List = PyTrilinos::pyDictToNewParameterList(obj);
    if (List == NULL) List = new Teuchos::ParameterList();
    // WARNING: THIS IS DELICATE, NULLSPACE SHOULD NOT DISAPPEAR
    // otherwise the pointer here stored will vanish. This function should
    // be used only through the MatrixPortal
    double* NullSpacePtr = (double*)NullSpace.Values();
    int NullSpaceDim = NullSpace.NumVectors();
    List->set("null space: type", "pre-computed");
    List->set("null space: vectors", NullSpacePtr);
    List->set("null space: dimension", NullSpaceDim);
    self->SetParameterList(*List);
    delete List;
    return(0);
  }
}
}
#endif

///////////////////////////////////
// MLAPI_InverseOperator support //
///////////////////////////////////
%include "MLAPI_InverseOperator.h"
namespace MLAPI
{
  %extend InverseOperator
  {
    MultiVector __mul__(const MultiVector& rhs) const
    {
      MultiVector result(rhs.GetVectorSpace(), rhs.GetNumVectors());
      self->Apply(rhs, result);
      return result;
    }
    bool Reshape(const Operator& Op, const std::string Type, PyObject* obj)
    {
      Teuchos::ParameterList * List = PyTrilinos::pyDictToNewParameterList(obj);
      if (List == NULL) return false;
      else
	self->Reshape(Op, Type, *List);
      delete List;
      return true;
    }
  }
}

//////////////////////////////////
// MLAPI_Operator_Utils support //
//////////////////////////////////
%include "MLAPI_Operator_Utils.h"

//////////////////////////////////////
// MLAPI_EpetraBaseOperator support //
//////////////////////////////////////
#ifdef HAVE_PYTRILINOS_EPETRA
%teuchos_rcp(MLAPI::EpetraBaseOperator)
%include "MLAPI_EpetraBaseOperator.h"
#endif

//////////////////////////
// MLAPI_Krylov support //
//////////////////////////
%include "MLAPI_Krylov.h"

///////////////////////////////
// MLAPI_Expressions support //
///////////////////////////////
%include "MLAPI_Expressions.h"

///////////////////////////////
// MLAPI_Aggregation support //
///////////////////////////////
%include "MLAPI_Aggregation.h"

///////////////////////
// MLAPI_Eig support //
///////////////////////
%include "MLAPI_Eig.h"

///////////////////////////
// MLAPI_Gallery support //
///////////////////////////
%include "MLAPI_Gallery.h"

%inline
{
  MLAPI::Operator GetPNonSmoothed(const MLAPI::Operator& A,
				  const MLAPI::MultiVector& ThisNS,
				  MLAPI::MultiVector& NextNS,
				  PyObject* obj)
  {
    Teuchos::ParameterList * List = PyTrilinos::pyDictToNewParameterList(obj);
    MLAPI::Operator Ptent;
    MLAPI::GetPtent(A, *List, ThisNS, Ptent, NextNS);
    delete List;
    return(Ptent);
  }

  bool Iterate(const MLAPI::Operator& A, const MLAPI::MultiVector& LHS,
	       const MLAPI::MultiVector& RHS, const MLAPI::BaseOperator& Prec, 
	       PyObject* obj)
  {
    Teuchos::ParameterList * List = PyTrilinos::pyDictToNewParameterList(obj);
    if (List == NULL) return(false);
    Krylov(A, LHS, RHS, Prec, *List);
    delete List;
    return(true);
  }
}

// Turn off the exception handling
%exception;
