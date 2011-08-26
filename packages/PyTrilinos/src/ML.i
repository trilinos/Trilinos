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
//  PyTrilinos includes
#include "PyTrilinos_PythonException.h"
#include "PyTrilinos_FILEstream.h"

// System includes
#include <iostream>
#include <sstream>
#include <vector>

// Configuration includes
#include "PyTrilinos_config.h"
#ifdef HAVE_INTTYPES_H
#undef HAVE_INTTYPES_H
#endif
#ifdef HAVE_STDINT_H
#undef HAVE_STDINT_H
#endif

// Teuchos includes
#ifdef HAVE_TEUCHOS
#include "Teuchos_RCPDecl.hpp"
#include "PyTrilinos_Teuchos_Util.h"
#endif

// Epetra includes
#ifdef HAVE_EPETRA
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_MapColoring.h"
#include "Epetra_FEVector.h"
#include "Epetra_Operator.h"
#include "Epetra_InvOperator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDistributor.h"
#include "Epetra_SerialDenseSVD.h"
#include "Epetra_Export.h"
#include "Epetra_OffsetIndex.h"

// Epetra python includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"
#include "Epetra_NumPyIntVector.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyFEVector.h"
#include "Epetra_NumPyIntSerialDenseMatrix.h"
#include "Epetra_NumPyIntSerialDenseVector.h"
#include "Epetra_NumPySerialDenseMatrix.h"
#include "Epetra_NumPySerialSymDenseMatrix.h"
#include "Epetra_NumPySerialDenseVector.h"
#endif

// IFPACK includes
#include "Ifpack_IC.h"
#include "Ifpack_ICT.h"
#include "Ifpack_ILU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_Amesos.h"

// ML includes
#undef HAVE_STDINT_H
#undef HAVE_INTTYPES_H
#undef HAVE_SYS_TIME_H
#include "ml_MultiLevelPreconditioner.h"
#include "MLAPI.h"
#include "PyTrilinos_ML_Util.h"

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
#ifdef HAVE_TEUCHOS
%import "Teuchos.i"
#endif
#ifdef HAVE_EPETRA
%include "Epetra_RowMatrix_Utils.i"
%ignore Epetra_Version;
%import  "Epetra.i"
%ignore Ifpack_Version;
%import  "IFPACK.i"
#endif

// Teuchos::RCP handling
#ifdef HAVE_TEUCHOS
%teuchos_rcp(MLAPI::DoubleVector)
%teuchos_rcp(MLAPI::ML_Operator_Box)
#ifdef HAVE_EPETRA
%teuchos_rcp(ML_Epetra::MultiLevelPreconditioner)
#endif
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
  SWIG_CATCH_STDEXCEPT
  catch(int ierr)
  {
    SWIG_exception(SWIG_RuntimeError, "ML exception thrown");
  }
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
      if (!PyFile_Check(ostream))
      {
	PyErr_SetString(PyExc_IOError, "Print() method expects file object");
	goto fail;
      }
      else
      {
	std::FILE * f = PyFile_AsFile(ostream);
	PyTrilinos::FILEstream buffer(f);
	std::ostream os(&buffer);
	self->Print(os, verbose);
	os.flush();
      }
    }
    return Py_BuildValue("");
  fail:
    return NULL;
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
%typemap(out) double* GetValues
{
  npy_intp dims = (arg1)->GetMyLength();
  $result = PyArray_SimpleNewFromData(1, &dims, PyArray_DOUBLE, (void*) $1);
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
      obj_to_array_contiguous_allow_conversion(pyValues, PyArray_DOUBLE,
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
#ifdef HAVE_EPETRA
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
#ifdef HAVE_EPETRA
%include "PyTrilinos_ML_Util.h"
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
#ifdef HAVE_TEUCHOS
#ifdef HAVE_EPETRA
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
#endif

///////////////////////////////////
// MLAPI_InverseOperator support //
///////////////////////////////////
%include "MLAPI_InverseOperator.h"
#ifdef HAVE_TEUCHOS
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
#endif

//////////////////////////////////
// MLAPI_Operator_Utils support //
//////////////////////////////////
%include "MLAPI_Operator_Utils.h"

//////////////////////////////////////
// MLAPI_EpetraBaseOperator support //
//////////////////////////////////////
#ifdef HAVE_EPETRA
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

#ifdef HAVE_TEUCHOS
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
#endif

// Turn off the exception handling
%exception;
