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
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseOperator.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialSymDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialDenseSVD.h"

%}

// This feature is turned on for most of the rest of the Epetra
// wrappers.  However, here it appears to cause us errors, so we turn
// it off.
%feature("compactdefaultargs", "0");

/////////////////////////////////////////////////////////
// Teuchos::RCP<> support for all classes in this file //
/////////////////////////////////////////////////////////
%teuchos_rcp(Epetra_IntSerialDenseMatrix)
%teuchos_rcp(Epetra_IntSerialDenseVector)
%teuchos_rcp(Epetra_SerialDenseOperator )
%teuchos_rcp(Epetra_SerialDenseMatrix   )
%teuchos_rcp(Epetra_SerialSymDenseMatrix)
%teuchos_rcp(Epetra_SerialDenseVector   )
%teuchos_rcp(Epetra_SerialDenseSolver   )
%teuchos_rcp(Epetra_SerialDenseSVD      )

/////////////////////////////////////////
// Epetra_IntSerialDenseMatrix support //
/////////////////////////////////////////
%rename(IntSerialDenseMatrix) Epetra_IntSerialDenseMatrix;
%inline
{
  PyObject *
    _extractNumPyArrayFromEpetraIntSerialDenseMatrix(
      const Epetra_IntSerialDenseMatrix & source)
  {
    // This NumPy function returns a borrowed pointer: do not DECREF
    PyArray_Descr * dtype = PyArray_DescrFromType(NPY_INT);
    npy_intp dim[2] = { source.M(), source.N() };
    int * data = const_cast< int* >(source.A());
    PyObject * result = PyArray_NewFromDescr(&PyArray_Type, dtype, 2, dim,
                                             NULL, (void*)data,
                                             NPY_ARRAY_FARRAY, NULL);
    return result;
  }
}
%extend Epetra_IntSerialDenseMatrix
{
  Epetra_IntSerialDenseMatrix(PyObject * array)
  {
    // This NumPy function returns a borrowed pointer: do not DECREF
    PyArray_Descr * dtype = PyArray_DescrFromType(NPY_INT);
    PyArrayObject * matrix =
      (PyArrayObject*)PyArray_FromAny(array, dtype, 2, 2,
                                      NPY_ARRAY_FARRAY, NULL);
    if (!matrix) throw PyTrilinos::PythonException();
    int nRows  = PyArray_DIM(matrix,0);
    int nCols  = PyArray_DIM(matrix,1);
    int * data = (int*)PyArray_DATA(matrix);
    Epetra_IntSerialDenseMatrix * result =
      new Epetra_IntSerialDenseMatrix(Copy, data, nRows, nRows, nCols);
    Py_DECREF(matrix);
    return result;
  }
}
%pythonappend Epetra_IntSerialDenseMatrix::Shape
{
    if "array" in self.__dict__: del self.__dict__["array"]
}
%pythonappend Epetra_IntSerialDenseMatrix::Reshape
{
    if "array" in self.__dict__: del self.__dict__["array"]
}
%ignore Epetra_IntSerialDenseMatrix::Epetra_IntSerialDenseMatrix(Epetra_DataAccess,
                                                                 int*,int,int,int);
%ignore Epetra_IntSerialDenseMatrix::operator()(int,int);
%ignore Epetra_IntSerialDenseMatrix::A() const;
%ignore Epetra_IntSerialDenseMatrix::MakeViewOf;
%include "Epetra_IntSerialDenseMatrix.h"
%pythoncode
%{
  def IntSerialDenseMatrix_getattr(self, name):
      if name == "array":
          a = _extractNumPyArrayFromEpetraIntSerialDenseMatrix(self)
          #a = a.transpose()
          self.__dict__["array"] = a
          return a
      elif name == "shape":
          return self.array.shape
      elif name == "dtype":
          return self.array.dtype
      raise AttributeError("'%s' not an attribute of IntSerialDenseMatrix" % name)
  def IntSerialDenseMatrix_setattr(self, name, value):
      if name == "array":
          raise AttributeError("Cannot change IntSerialDenseMatrix 'array' attribute")
      elif name == "shape":
          self.array.shape = value
      elif name == "dtype":
          raise AttributeError("Cannot change IntSerialDenseMatrix 'dtype' attribute")
      else:
          self.__dict__[name] = value
  IntSerialDenseMatrix.__getattr__ = IntSerialDenseMatrix_getattr
  IntSerialDenseMatrix.__setattr__ = IntSerialDenseMatrix_setattr
  IntSerialDenseMatrix.__getitem__ = lambda self, i: self.array.__getitem__(i)
  IntSerialDenseMatrix.__setitem__ = lambda self, i, v: self.array.__setitem__(i,v)
  IntSerialDenseMatrix.__len__     = lambda self: self.array.__len__()
  IntSerialDenseMatrix.__str__     = lambda self: self.array.__str__()
  IntSerialDenseMatrix.A           = lambda self: self.array
  class_array_add_math(IntSerialDenseMatrix)
  class_array_add_comp(IntSerialDenseMatrix)
%}

/////////////////////////////////////////
// Epetra_IntSerialDenseVector support //
/////////////////////////////////////////
%rename(IntSerialDenseVector) Epetra_IntSerialDenseVector;
%inline
{
  PyObject *
    _extractNumPyArrayFromEpetraIntSerialDenseVector(
      const Epetra_IntSerialDenseVector & source)
  {
    npy_intp dim[1] = { source.Length() };
    int * data = const_cast< int* >(source.Values());
    return PyArray_SimpleNewFromData(1, dim, NPY_INT, (void*)data);
  }
}
%extend Epetra_IntSerialDenseVector
{
  Epetra_IntSerialDenseVector(PyObject * array)
  {
    PyArrayObject * vector =
      (PyArrayObject*) PyArray_ContiguousFromAny(array, NPY_INT, 0, 0);
    if (!vector) throw PyTrilinos::PythonException();
    int ndim   = PyArray_NDIM(vector);
    int size = 1;
    for (int i=0; i<ndim; ++i) size *= PyArray_DIM(vector,i);
    int * data = (int*)PyArray_DATA(vector);
    Epetra_IntSerialDenseVector * result =
      new Epetra_IntSerialDenseVector(Copy, data, size);
    Py_DECREF(vector);
    return result;
  }
}
%pythonappend Epetra_IntSerialDenseVector::Size
{
    if "array" in self.__dict__: del self.__dict__["array"]
}
%pythonappend Epetra_IntSerialDenseVector::Resize
{
    if "array" in self.__dict__: del self.__dict__["array"]
}
%ignore Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector(Epetra_DataAccess,
                                                                 int*,int);
%ignore Epetra_IntSerialDenseVector::operator()(int);
%ignore Epetra_IntSerialDenseVector::Values;
%include "Epetra_IntSerialDenseVector.h"
%pythoncode
{
  def IntSerialDenseVector_getattr(self, name):
      if name == "array":
          a = _extractNumPyArrayFromEpetraIntSerialDenseVector(self)
          self.__dict__["array"] = a
          return a
      elif name == "shape":
          return self.array.shape
      elif name == "dtype":
          return self.array.dtype
      raise AttributeError("'%s' not an attribute of IntSerialDenseVector" % name)
  def IntSerialDenseVector_setattr(self, name, value):
      if name == "array":
          raise AttributeError("Cannot change IntSerialDenseVector 'array' attribute")
      elif name == "shape":
          self.array.shape = value
      elif name == "dtype":
          raise AttributeError("Cannot change IntSerialDenseVector 'dtype' attribute")
      else:
          self.__dict__[name] = value
  IntSerialDenseVector.__getattr__ = IntSerialDenseVector_getattr
  IntSerialDenseVector.__setattr__ = IntSerialDenseVector_setattr
  IntSerialDenseVector.__getitem__ = lambda self, i: self.array.__getitem__(i)
  IntSerialDenseVector.__setitem__ = lambda self, i, v: self.array.__setitem__(i,v)
  IntSerialDenseVector.__len__     = lambda self: self.array.__len__()
  IntSerialDenseVector.__str__     = lambda self: self.array.__str__()
  IntSerialDenseVector.Values      = lambda self: self.array
  class_array_add_math(IntSerialDenseVector)
  class_array_add_comp(IntSerialDenseVector)

}

////////////////////////////////////////
// Epetra_SerialDenseOperator support //
////////////////////////////////////////
%rename(SerialDenseOperator) Epetra_SerialDenseOperator;
%include "Epetra_SerialDenseOperator.h"

//////////////////////////////////////
// Epetra_SerialDenseMatrix support //
//////////////////////////////////////
%rename(SerialDenseMatrix) Epetra_SerialDenseMatrix;
%inline
{
  PyObject *
    _extractNumPyArrayFromEpetraSerialDenseMatrix(
      const Epetra_SerialDenseMatrix & source)
  {
    // This NumPy function returns a borrowed pointer: do not DECREF
    PyArray_Descr * dtype = PyArray_DescrFromType(NPY_DOUBLE);
    npy_intp dim[2] = { source.M(), source.N() };
    double * data = const_cast< double* >(source.A());
    PyObject * result = PyArray_NewFromDescr(&PyArray_Type, dtype, 2, dim,
                                             NULL, (void*)data,
                                             NPY_ARRAY_FARRAY, NULL);
    return result;
  }
}
%extend Epetra_SerialDenseMatrix
{
  Epetra_SerialDenseMatrix(PyObject * array,
                           bool set_object_label = true)
  {
    // This NumPy function returns a borrowed pointer: do not DECREF
    PyArray_Descr * dtype = PyArray_DescrFromType(NPY_DOUBLE);
    PyArrayObject * matrix =
      (PyArrayObject*) PyArray_FromAny(array, dtype, 2, 2,
                                       NPY_ARRAY_FARRAY, NULL);
    if (!matrix) throw PyTrilinos::PythonException();
    int nRows     = PyArray_DIM(matrix,0);
    int nCols     = PyArray_DIM(matrix,1);
    double * data = (double*)PyArray_DATA(matrix);
    Epetra_SerialDenseMatrix * result =
      new Epetra_SerialDenseMatrix(Copy, data, nRows, nRows, nCols,
                                   set_object_label);
    Py_DECREF(matrix);
    return result;
  }
}
%pythonappend Epetra_SerialDenseMatrix::Shape
{
    if "array" in self.__dict__: del self.__dict__["array"]
}
%pythonappend Epetra_SerialDenseMatrix::Reshape
{
    if "array" in self.__dict__: del self.__dict__["array"]
}
%ignore Epetra_SerialDenseMatrix::operator()(int,int);
%ignore Epetra_SerialDenseMatrix::A() const;
%include "Epetra_SerialDenseMatrix.h"
%pythoncode
%{
  def SerialDenseMatrix_getattr(self, name):
      if name == "array":
          a = _extractNumPyArrayFromEpetraSerialDenseMatrix(self)
          #a = a.transpose()
          self.__dict__["array"] = a
          return a
      elif name == "shape":
          return self.array.shape
      elif name == "dtype":
          return self.array.dtype
      raise AttributeError("'%s' not an attribute of SerialDenseMatrix" % name)
  def SerialDenseMatrix_setattr(self, name, value):
      if name == "array":
          raise AttributeError("Cannot change SerialDenseMatrix 'array' attribute")
      elif name == "shape":
          self.array.shape = value
      elif name == "dtype":
          raise AttributeError("Cannot change SerialDenseMatrix 'dtype' attribute")
      else:
          self.__dict__[name] = value
  SerialDenseMatrix.__getattr__ = SerialDenseMatrix_getattr
  SerialDenseMatrix.__setattr__ = SerialDenseMatrix_setattr
  SerialDenseMatrix.__getitem__ = lambda self, i: self.array.__getitem__(i)
  SerialDenseMatrix.__setitem__ = lambda self, i, v: self.array.__setitem__(i,v)
  SerialDenseMatrix.__len__     = lambda self: self.array.__len__()
  SerialDenseMatrix.__str__     = lambda self: self.array.__str__()
  SerialDenseMatrix.A           = lambda self: self.array
  class_array_add_math(SerialDenseMatrix)
  class_array_add_comp(SerialDenseMatrix)

%}

/////////////////////////////////////////
// Epetra_SerialSymDenseMatrix support //
/////////////////////////////////////////
%rename(SerialSymDenseMatrix) Epetra_SerialSymDenseMatrix;
%include "Epetra_SerialSymDenseMatrix.h"

//////////////////////////////////////
// Epetra_SerialDenseVector support //
//////////////////////////////////////
%rename(SerialDenseVector) Epetra_SerialDenseVector;
%inline
{
  PyObject *
    _extractNumPyArrayFromEpetraSerialDenseVector(
      const Epetra_SerialDenseVector & source)
  {
    npy_intp dim[1] = { source.Length() };
    double * data = const_cast< double* >(source.Values());
    return PyArray_SimpleNewFromData(1, dim, NPY_DOUBLE, (void*)data);
  }
}
%extend Epetra_SerialDenseVector
{
  Epetra_SerialDenseVector(PyObject * array)
  {
    PyArrayObject * vector =
      (PyArrayObject*) PyArray_ContiguousFromAny(array, NPY_DOUBLE, 0, 0);
    if (!vector) throw PyTrilinos::PythonException();
    int ndim   = PyArray_NDIM(vector);
    int size = 1;
    for (int i=0; i<ndim; ++i) size *= PyArray_DIM(vector,i);
    double * data = (double*)PyArray_DATA(vector);
    Epetra_SerialDenseVector * result =
      new Epetra_SerialDenseVector(Copy, data, size);
    Py_DECREF(vector);
    return result;
  }
}
%pythonappend Epetra_SerialDenseVector::Size
{
    if "array" in self.__dict__: del self.__dict__["array"]
}
%pythonappend Epetra_SerialDenseVector::Resize
{
    if "array" in self.__dict__: del self.__dict__["array"]
}
%ignore Epetra_SerialDenseVector::Epetra_SerialDenseVector(Epetra_DataAccess,
                                                           double*,int);
%ignore Epetra_SerialDenseVector::operator()(int);
%include "Epetra_SerialDenseVector.h"
%pythoncode
%{
  def SerialDenseVector_getattr(self, name):
      if name == "array":
          a = _extractNumPyArrayFromEpetraSerialDenseVector(self)
          self.__dict__["array"] = a
          return a
      elif name == "shape":
          return self.array.shape
      elif name == "dtype":
          return self.array.dtype
      raise AttributeError("'%s' not an attribute of SerialDenseVector" % name)
  def SerialDenseVector_setattr(self, name, value):
      if name == "array":
          raise AttributeError("Cannot change SerialDenseVector 'array' attribute")
      elif name == "shape":
          self.array.shape = value
      elif name == "dtype":
          raise AttributeError("Cannot change SerialDenseVector 'dtype' attribute")
      else:
          self.__dict__[name] = value
  SerialDenseVector.__getattr__ = SerialDenseVector_getattr
  SerialDenseVector.__setattr__ = SerialDenseVector_setattr
  SerialDenseVector.__getitem__ = lambda self, i: self.array.__getitem__(i)
  SerialDenseVector.__setitem__ = lambda self, i, v: self.array.__setitem__(i,v)
  SerialDenseVector.__len__     = lambda self: self.array.__len__()
  SerialDenseVector.__str__     = lambda self: self.array.__str__()
  SerialDenseVector.Values      = lambda self: self.array
  class_array_add_math(SerialDenseVector)
  class_array_add_comp(SerialDenseVector)

%}

//////////////////////////////////////
// Epetra_SerialDenseSolver support //
//////////////////////////////////////
%ignore Epetra_SerialDenseSolver::ReciprocalConditionEstimate(double&);
%rename(SerialDenseSolver) Epetra_SerialDenseSolver;
%fragment("NumPy_Macros");  // These macros depend upon this fragment
%epetra_intarray1d_output_method(Epetra_SerialDenseSolver,IPIV,M)
%epetra_array2d_output_method(Epetra_SerialDenseSolver,A,M,N    )
%epetra_array2d_output_method(Epetra_SerialDenseSolver,B,N,NRHS )
%epetra_array2d_output_method(Epetra_SerialDenseSolver,X,N,NRHS )
%epetra_array2d_output_method(Epetra_SerialDenseSolver,AF,M,N   )
%epetra_array1d_output_method(Epetra_SerialDenseSolver,FERR,NRHS)
%epetra_array1d_output_method(Epetra_SerialDenseSolver,BERR,NRHS)
%epetra_array1d_output_method(Epetra_SerialDenseSolver,R,M      )
%epetra_array1d_output_method(Epetra_SerialDenseSolver,C,N      )
%extend Epetra_SerialDenseSolver
{
  double ReciprocalConditionEstimate()
  {
    double value = 0.0;
    int result = self->ReciprocalConditionEstimate(value);
    if (result)
    {
      PyErr_Format(PyExc_RuntimeError,
		   "ReciprocalConditionEstimate method returned LAPACK error code %d",
		   result);
      value = -1.0;
    }
    return value;
  }
}
%include "Epetra_SerialDenseSolver.h"

///////////////////////////////////
// Epetra_SerialDenseSVD support //
///////////////////////////////////
%rename(SerialDenseSVD) Epetra_SerialDenseSVD;
%include "Epetra_SerialDenseSVD.h"

%feature("compactdefaultargs");       // Turn the feature back on
