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
#include "PyTrilinos_PythonException.hpp"

// Epetra include files
#include "Epetra_SerialComm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntVector.h"

%}

// PyTrilinos configuration
%include "PyTrilinos_config.h"

/////////////////////////////////////////////////////////
// Teuchos::RCP<> support for all classes in this file //
/////////////////////////////////////////////////////////
%teuchos_rcp_dap(PyTrilinos::convertPythonToEpetraIntVector  , Epetra_IntVector  )
%teuchos_rcp_dap(PyTrilinos::convertPythonToEpetraMultiVector, Epetra_MultiVector)
%teuchos_rcp_dap(PyTrilinos::convertPythonToEpetraVector     , Epetra_Vector     )
%teuchos_rcp(Epetra_FEVector   )
%teuchos_rcp_epetra_argout(Epetra_IntVector  )
%teuchos_rcp_epetra_argout(Epetra_MultiVector)
%teuchos_rcp_epetra_argout(Epetra_Vector     )
%teuchos_rcp_epetra_argout(Epetra_FEVector   )

//////////////////////////////
// Epetra_IntVector support //
//////////////////////////////
%rename(IntVector) Epetra_IntVector;
%inline
{
  PyObject *
    _extractNumPyArrayFromEpetraIntVector(const Epetra_IntVector & source)
  {
    npy_intp dim[1] = { source.MyLength() };
    int * data = source.Values();
    return PyArray_SimpleNewFromData(1, dim, NPY_INT, (void*)data);
  }
}
%extend Epetra_IntVector
{
  Epetra_IntVector(Epetra_DataAccess cv,
                   Epetra_BlockMap const & map,
                   PyObject * values)
  {
    PyObject * array = PyArray_ContiguousFromAny(values, NPY_INT, 0, 0);
    if (!array) throw PyTrilinos::PythonException();
    Epetra_IntVector * result =
      new Epetra_IntVector(cv,
                           map,
                           (int*)PyArray_DATA((PyArrayObject*)array));
    Py_DECREF(array);
    return result;
  }

  Epetra_IntVector(Epetra_BlockMap const & map,
                   PyObject * values)
  {
    PyObject * array = PyArray_ContiguousFromAny(values, NPY_INT, 0, 0);
    if (!array) throw PyTrilinos::PythonException();
    Epetra_IntVector * result =
      new Epetra_IntVector(Copy,
                           map,
                           (int*)PyArray_DATA((PyArrayObject*)array));
    Py_DECREF(array);
    return result;
  }

  Epetra_IntVector(PyObject * values)
  {
    Epetra_SerialComm comm;
    PyObject * array = PyArray_ContiguousFromAny(values, NPY_INT, 0, 0);
    if (!array) throw PyTrilinos::PythonException();
    int numPoints = 1;
    for (int i=0; i < PyArray_NDIM((PyArrayObject*)array); ++i)
      numPoints *= PyArray_DIM((PyArrayObject*)array,i);
    Epetra_Map map(numPoints, 0, comm);
    Epetra_IntVector * result =
      new Epetra_IntVector(Copy,
                           map,
                           (int*)PyArray_DATA((PyArrayObject*)array));
    Py_DECREF(array);
    return result;
  }

  PyObject * __distarray__()
  {
    return PyTrilinos::convertToDistArray(*self);
  }
}
%ignore Epetra_IntVector::Values;
%ignore Epetra_IntVector::Epetra_IntVector(Epetra_DataAccess,
                                           Epetra_BlockMap const &,
                                           int*);
%include "Epetra_IntVector.h"
%pythoncode
%{
  def IntVector_getattr(self, name):
      if name == "array":
          a = _extractNumPyArrayFromEpetraIntVector(self)
          self.__dict__["array"] = a
          return a
      elif name == "shape":
          return self.array.shape
      elif name == "dtype":
          return self.array.dtype
      raise AttributeError("'%s' not an attribute of IntVector" % name)
  def IntVector_setattr(self, name, value):
      if name == "array":
          raise AttributeError("Cannot change IntVector 'array' attribute")
      elif name == "shape":
          self.array.shape = value
      elif name == "dtype":
          raise AttributeError("Cannot change IntVector 'dtype' attribute")
      else:
          self.__dict__[name] = value
  IntVector.__getattr__ = IntVector_getattr
  IntVector.__setattr__ = IntVector_setattr
  IntVector.__getitem__ = lambda self, i: self.array.__getitem__(i)
  IntVector.__setitem__ = lambda self, i, v: self.array.__setitem__(i,v)
  IntVector.__len__     = lambda self: self.array.__len__()
  IntVector.__str__     = lambda self: self.array.__str__()
  IntVector.copy        = lambda self: IntVector(self)
  IntVector.ExtractCopy = lambda self: self.array.copy()
  IntVector.ExtractView = lambda self: self.array
  IntVector.Values      = lambda self: self.array
  class_array_add_math(IntVector)
  class_array_add_comp(IntVector)
%}

////////////////////////////////
// Epetra_MultiVector support //
////////////////////////////////
%rename(MultiVector) Epetra_MultiVector;
%inline
{
  PyObject *
    _extractNumPyArrayFromEpetraMultiVector(const Epetra_MultiVector & source)
  {
    npy_intp dims[2] = { source.NumVectors(), source.MyLength() };
    double * data = source.Values();
    PyObject * array = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, (void*)data);
    return array;
  }
}
%extend Epetra_MultiVector
{
  Epetra_MultiVector(Epetra_DataAccess cv,
                     const Epetra_BlockMap & map,
                     PyObject * values,
                     int numVectors)
  {
    PyObject * array = PyArray_ContiguousFromAny(values, NPY_DOUBLE, 2, 2);
    if (!array) throw PyTrilinos::PythonException();
    int arrayDim1 = PyArray_DIM((PyArrayObject*)array, 0);
    int arrayDim2 = PyArray_DIM((PyArrayObject*)array, 1);
    if (arrayDim1 != numVectors)
    {
      PyErr_Format(PyExc_ValueError,
                   "First array dimension of %d not equal to number of "
                   "vectors %d",
                   arrayDim1, numVectors);
      Py_DECREF(array);
      throw PyTrilinos::PythonException();
    }
    Epetra_MultiVector * result =
      new Epetra_MultiVector(cv,
                             map,
                             (double*)PyArray_DATA((PyArrayObject*)array),
                             arrayDim2,
                             arrayDim1);
    Py_DECREF(array);
    return result;
  }

  Epetra_MultiVector(const Epetra_BlockMap & map,
                     PyObject * obj)
  {
    PyArrayObject * array =
      (PyArrayObject*) PyArray_ContiguousFromObject(obj, NPY_DOUBLE, 0, 0);
    if (!array)
    {
      PyObject * type   = PyObject_Type(obj);
%#if PY_VERSION_HEX >= 0x03000000
      PyObject * string = PyUnicode_AsASCIIString(PyObject_Str(type));
%#else
      PyObject * string = PyObject_Str(type);
%#endif
      PyErr_Format(PyExc_TypeError,
                   "Epetra_MultiVector constructor expected an array-like "
                   "object as second argument, got %s",
                   PyString_AsString(string));
      Py_DECREF(type);
      Py_DECREF(string);
      throw PyTrilinos::PythonException();
    }
    int numVectors, numPoints;
    int ndim = PyArray_NDIM(array);
    if (ndim == 1)
    {
      numVectors = 1;
      numPoints  = PyArray_DIM(array, 0);
    }
    else
    {
      numVectors = PyArray_DIM(array, 0);
      numPoints  = 1;
      for (int i=1; i < ndim; ++i) numPoints *= PyArray_DIM(array, i);
    }
    if (numPoints < map.NumMyElements())
    {
      PyErr_Format(PyExc_ValueError,
                   "Epetra_MultiVector constructor expected second argument "
                   "to be an array with at least %d elements per vector, got "
                   "%d",
                   map.NumMyElements(), numPoints);
      throw PyTrilinos::PythonException();
    }
    if (!PyArray_ISFLOAT(array))
    {
      PyErr_SetString(PyExc_TypeError,
                      "Epetra_MultiVector constructor expected a NumPy array "
                      "of type double");
      throw PyTrilinos::PythonException();
    }
    double * data  = (double*)PyArray_DATA(array);
    Epetra_MultiVector * result =
      new Epetra_MultiVector(Copy, map, data, numPoints, numVectors);
    return result;
  }

  Epetra_MultiVector(PyObject * obj)
  {
    PyArrayObject * array =
      (PyArrayObject*) PyArray_ContiguousFromObject(obj, NPY_DOUBLE, 0, 0);
    if (!array)
    {
      PyObject * type   = PyObject_Type(obj);
%#if PY_VERSION_HEX >= 0x03000000
      PyObject * string = PyUnicode_AsASCIIString(PyObject_Str(type));
%#else
      PyObject * string = PyObject_Str(type);
%#endif
      PyErr_Format(PyExc_TypeError,
                   "Epetra_MultiVector constructor expected an array-like "
                   "object as first argument, got %s",
                   PyString_AsString(string));
      Py_DECREF(type);
      Py_DECREF(string);
      throw PyTrilinos::PythonException();
    }
    int numVectors, numPoints;
    int ndim = PyArray_NDIM(array);
    if (ndim == 1)
    {
      numVectors = 1;
      numPoints  = PyArray_DIM(array, 0);
    }
    else
    {
      numVectors = PyArray_DIM(array, 0);
      numPoints  = 1;
      for (int i=1; i < ndim; ++i) numPoints *= PyArray_DIM(array, i);
    }
    Epetra_SerialComm comm;
    Epetra_Map map(numPoints,0,comm);
    if (!PyArray_ISFLOAT(array))
    {
      PyErr_SetString(PyExc_TypeError,
                      "Epetra_MultiVector constructor expected a NumPy array "
                      "of type double");
      throw PyTrilinos::PythonException();
    }
    double * data  = (double*)PyArray_DATA(array);
    Epetra_MultiVector * result =
      new Epetra_MultiVector(Copy, map, data, numPoints, numVectors);
    return result;
  }

  Epetra_MultiVector(Epetra_DataAccess cv,
                     const Epetra_MultiVector & source)
  {
    return new Epetra_MultiVector(cv, source, 0, source.NumVectors());
  }

  PyObject * Dot(const Epetra_MultiVector & a) const
  {
    npy_intp dims[1] = { self->NumVectors() };
    PyObject * result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    double * data = (double*)PyArray_DATA((PyArrayObject*)result);
    self->Dot(a, data);
    return result;
  }

  PyObject * Norm1() const
  {
    npy_intp dims[1] = { self->NumVectors() };
    PyObject * result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    double * data = (double*)PyArray_DATA((PyArrayObject*)result);
    self->Norm1(data);
    return result;
  }

  PyObject * Norm2() const
  {
    npy_intp dims[1] = { self->NumVectors() };
    PyObject * result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    double * data = (double*)PyArray_DATA((PyArrayObject*)result);
    self->Norm2(data);
    return result;
  }

  PyObject * NormInf() const
  {
    npy_intp dims[1] = { self->NumVectors() };
    PyObject * result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    double * data = (double*)PyArray_DATA((PyArrayObject*)result);
    self->NormInf(data);
    return result;
  }

  PyObject * NormWeighted(const Epetra_MultiVector & weights) const
  {
    npy_intp dims[1] = { self->NumVectors() };
    PyObject * result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    double * data = (double*)PyArray_DATA((PyArrayObject*)result);
    self->NormWeighted(weights, data);
    return result;
  }

  PyObject * MinValue() const
  {
    npy_intp dims[1] = { self->NumVectors() };
    PyObject * result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    double * data = (double*)PyArray_DATA((PyArrayObject*)result);
    self->MinValue(data);
    return result;
  }

  PyObject * MaxValue() const
  {
    npy_intp dims[1] = { self->NumVectors() };
    PyObject * result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    double * data = (double*)PyArray_DATA((PyArrayObject*)result);
    self->MaxValue(data);
    return result;
  }

  PyObject * MeanValue() const
  {
    npy_intp dims[1] = { self->NumVectors() };
    PyObject * result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    double * data = (double*)PyArray_DATA((PyArrayObject*)result);
    self->MeanValue(data);
    return result;
  }

  const Epetra_Vector & __call__(int i) const
  {
    const Epetra_Vector & result = *(self->operator()(i));
    return result;
  }

  PyObject * __distarray__() const
  {
    return PyTrilinos::convertToDistArray(*self);
  }
}
%ignore Epetra_MultiVector::Epetra_MultiVector(Epetra_DataAccess,
                                               const Epetra_BlockMap &,
                                               double *,
                                               int,
                                               int);
%ignore Epetra_MultiVector::Epetra_MultiVector(Epetra_DataAccess,
                                               const Epetra_BlockMap &,
                                               double **,
                                               int,
                                               int);
%ignore Epetra_MultiVector::operator()(int);
%ignore Epetra_MultiVector::operator()(int) const;
%ignore Epetra_MultiVector::ExtractCopy(double *, int   ) const;
%ignore Epetra_MultiVector::ExtractCopy(double **       ) const;
%ignore Epetra_MultiVector::ExtractView(double **, int *) const;
%ignore Epetra_MultiVector::ExtractView(double ***      ) const;
%ignore Epetra_MultiVector::Dot(const Epetra_MultiVector&,double*) const;
%ignore Epetra_MultiVector::Norm1(double*) const;
%ignore Epetra_MultiVector::Norm2(double*) const;
%ignore Epetra_MultiVector::NormInf(double*) const;
%ignore Epetra_MultiVector::NormWeighted(const Epetra_MultiVector&,double*) const;
%ignore Epetra_MultiVector::MinValue(double*) const;
%ignore Epetra_MultiVector::MaxValue(double*) const;
%ignore Epetra_MultiVector::MeanValue(double*) const;
%ignore Epetra_MultiVector::ResetView(double **);
%ignore Epetra_MultiVector::Pointers() const;
%apply (int* IN_ARRAY1, int DIM1) {(int* Indices, int NumVectors)};
%include "Epetra_MultiVector.h"
%pythoncode
%{
  def MultiVector_getattr(self, name):
      if name == "array":
          a = _extractNumPyArrayFromEpetraMultiVector(self)
          self.__dict__["array"] = a
          return a
      elif name == "shape":
          return self.array.shape
      elif name == "dtype":
          return self.array.dtype
      raise AttributeError("'%s' not an attribute of MultiVector" % name)
  def MultiVector_setattr(self, name, value):
      if name == "array":
          raise AttributeError("Cannot change MultiVector 'array' attribute")
      elif name == "shape":
          self.array.shape = value
      elif name == "dtype":
          raise AttributeError("Cannot change MultiVector 'dtype' attribute")
      else:
          self.__dict__[name] = value
  def MultiVector_getitem(self,i):
      if isinstance(i,int):
          return self.__call__(i)
      else:
          return self.array.__getitem__(i)
  MultiVector.__getattr__ = MultiVector_getattr
  MultiVector.__setattr__ = MultiVector_setattr
  MultiVector.__getitem__ = MultiVector_getitem
  MultiVector.__setitem__ = lambda self, i, v: self.array.__setitem__(i,v)
  MultiVector.__len__     = lambda self: self.array.__len__()
  MultiVector.__str__     = lambda self: self.array.__str__()
  MultiVector.copy        = lambda self: MultiVector(self)
  MultiVector.ExtractCopy = lambda self: self.array.copy()
  MultiVector.ExtractView = lambda self: self.array
  class_array_add_math(MultiVector)
  class_array_add_comp(MultiVector)
%}

///////////////////////////
// Epetra_Vector support //
///////////////////////////
%rename(Vector) Epetra_Vector;
%inline
{
  PyObject *
    _extractNumPyArrayFromEpetraVector(const Epetra_Vector & source)
  {
    npy_intp dim[1] = { source.MyLength() };
    double * data = source.Values();
    return PyArray_SimpleNewFromData(1, dim, NPY_DOUBLE, (void*)data);
  }
}
%extend Epetra_Vector
{
  Epetra_Vector(Epetra_DataAccess cv,
                Epetra_BlockMap const & map,
                PyObject * values)
  {
    Epetra_Vector * result = 0;
    PyObject * array = PyArray_ContiguousFromAny(values, NPY_DOUBLE, 1, 1);
    if (!array) SWIG_fail;
    result = new Epetra_Vector(cv,
                               map,
                               (double*)PyArray_DATA((PyArrayObject*)array));
    Py_DECREF(array);
    return result;
  fail:
    Py_XDECREF(array);
    return NULL;
  }

  Epetra_Vector(const Epetra_BlockMap & map,
                PyObject * obj)
  {
    PyArrayObject * array = NULL;
    if (PyArray_Check(obj))
    {
      array = (PyArrayObject*)obj;
    }
    else
    {
      array =
        (PyArrayObject*) PyArray_ContiguousFromObject(obj, NPY_DOUBLE, 0, 0);
      if (!array)
      {
        PyObject * type   = PyObject_Type(obj);
%#if PY_VERSION_HEX >= 0x03000000
        PyObject * string = PyUnicode_AsASCIIString(PyObject_Str(type));
%#else
        PyObject * string = PyObject_Str(type);
%#endif
        PyErr_Format(PyExc_TypeError,
                     "Epetra_Vector constructor expected an array-like object "
                     "as second argument, got %s",
                     PyString_AsString(string));
        Py_DECREF(type);
        Py_DECREF(string);
        throw PyTrilinos::PythonException();
      }
    }
    int ndim      = PyArray_NDIM(array);
    int numPoints = 1;
    for (int i=0; i < ndim; ++i) numPoints *= PyArray_DIM(array, i);
    if (numPoints < map.NumMyElements())
    {
      PyErr_Format(PyExc_ValueError,
                   "Epetra_Vector constructor expected second argument to be "
                   "an array with at least %d elements, got %d",
                   map.NumMyElements(), numPoints);
      throw PyTrilinos::PythonException();
    }
    if (!PyArray_ISFLOAT(array))
    {
      PyErr_SetString(PyExc_TypeError,
                      "Epetra_Vector constructor expected a NumPy array of "
                      "type double");
      throw PyTrilinos::PythonException();
    }
    double * data  = (double*)PyArray_DATA(array);
    Epetra_Vector * result = new Epetra_Vector(Copy, map, data);
    return result;
  }

  Epetra_Vector(PyObject * obj)
  {
    PyArrayObject * array =
      (PyArrayObject*) PyArray_ContiguousFromObject(obj, NPY_DOUBLE, 0, 0);
    if (!array)
    {
      PyObject * type   = PyObject_Type(obj);
%#if PY_VERSION_HEX >= 0x03000000
      PyObject * string = PyUnicode_AsASCIIString(PyObject_Str(type));
%#else
      PyObject * string = PyObject_Str(type);
%#endif
      PyErr_Format(PyExc_TypeError,
                   "Epetra_Vector constructor expected an array-like "
                   "object as first argument, got %s",
                   PyString_AsString(string));
      Py_DECREF(type);
      Py_DECREF(string);
      throw PyTrilinos::PythonException();
    }
    int ndim      = PyArray_NDIM(array);
    int numPoints = 1;
    for (int i=0; i < ndim; ++i) numPoints *= PyArray_DIM(array, i);
    Epetra_SerialComm comm;
    Epetra_Map map(numPoints,0,comm);
    if (!PyArray_ISFLOAT(array))
    {
      PyErr_SetString(PyExc_TypeError,
                      "Epetra_Vector constructor expected a NumPy array "
                      "of type double");
      throw PyTrilinos::PythonException();
    }
    double * data = (double*)PyArray_DATA(array);
    Epetra_Vector * result = new Epetra_Vector(Copy, map, data);
    return result;
  }

  Epetra_Vector(Epetra_DataAccess cv,
                const Epetra_Vector & source)
  {
    return new Epetra_Vector(cv, source, 0);
  }

  int ReplaceGlobalValues(PyObject * values,
                          PyObject * indices)
  {
    PyObject * vArray = PyArray_ContiguousFromAny(values, NPY_DOUBLE, 1, 1);
    if (!vArray) throw PyTrilinos::PythonException();
    PyObject * iArray = PyArray_ContiguousFromAny(indices, NPY_INT, 1, 1);
    if (!iArray)
    {
      Py_DECREF(vArray);
      throw PyTrilinos::PythonException();
    }
    int vLen = PyArray_DIM((PyArrayObject*)vArray, 0);
    int iLen = PyArray_DIM((PyArrayObject*)iArray, 0);
    if (vLen != iLen)
    {
      Py_DECREF(vArray);
      Py_DECREF(iArray);
      PyErr_Format(PyExc_ValueError,
                   "Length of values array (%d) does not equal length of "
                   "indices array (%d)",
                   vLen, iLen);
      throw PyTrilinos::PythonException();
    }
    double * vData = (double*) PyArray_DATA((PyArrayObject*)vArray);
    int *    iData = (int*)    PyArray_DATA((PyArrayObject*)iArray);
    int result = self->ReplaceGlobalValues(vLen, vData, iData);
    Py_DECREF(vArray);
    Py_DECREF(iArray);
    return result;
  }

  int ReplaceGlobalValues(int blockOffset,
                          PyObject * values,
                          PyObject * indices)
  {
    PyObject * vArray = PyArray_ContiguousFromAny(values, NPY_DOUBLE, 1, 1);
    if (!vArray) throw PyTrilinos::PythonException();
    PyObject * iArray = PyArray_ContiguousFromAny(indices, NPY_INT, 1, 1);
    if (!iArray)
    {
      Py_DECREF(vArray);
      throw PyTrilinos::PythonException();
    }
    int vLen = PyArray_DIM((PyArrayObject*)vArray, 0);
    int iLen = PyArray_DIM((PyArrayObject*)iArray, 0);
    if (vLen != iLen)
    {
      Py_DECREF(vArray);
      Py_DECREF(iArray);
      PyErr_Format(PyExc_ValueError,
                   "Length of values array (%d) does not equal length of "
                   "indices array (%d)",
                   vLen, iLen);
      throw PyTrilinos::PythonException();
    }
    double * vData = (double*) PyArray_DATA((PyArrayObject*)vArray);
    int *    iData = (int*)    PyArray_DATA((PyArrayObject*)iArray);
    int result = self->ReplaceGlobalValues(vLen, blockOffset, vData, iData);
    Py_DECREF(vArray);
    Py_DECREF(iArray);
    return result;
  }

  int ReplaceMyValues(PyObject * values,
                      PyObject * indices)
  {
    PyObject * vArray = PyArray_ContiguousFromAny(values, NPY_DOUBLE, 1, 1);
    if (!vArray) throw PyTrilinos::PythonException();
    PyObject * iArray = PyArray_ContiguousFromAny(indices, NPY_INT, 1, 1);
    if (!iArray)
    {
      Py_DECREF(vArray);
      throw PyTrilinos::PythonException();
    }
    int vLen = PyArray_DIM((PyArrayObject*)vArray, 0);
    int iLen = PyArray_DIM((PyArrayObject*)iArray, 0);
    if (vLen != iLen)
    {
      Py_DECREF(vArray);
      Py_DECREF(iArray);
      PyErr_Format(PyExc_ValueError,
                   "Length of values array (%d) does not equal length of "
                   "indices array (%d)",
                   vLen, iLen);
      throw PyTrilinos::PythonException();
    }
    double * vData = (double*) PyArray_DATA((PyArrayObject*)vArray);
    int *    iData = (int*)    PyArray_DATA((PyArrayObject*)iArray);
    int result = self->ReplaceMyValues(vLen, vData, iData);
    Py_DECREF(vArray);
    Py_DECREF(iArray);
    return result;
  }

  int ReplaceMyValues(int blockOffset,
                      PyObject * values,
                      PyObject * indices)
  {
    PyObject * vArray = PyArray_ContiguousFromAny(values, NPY_DOUBLE, 1, 1);
    if (!vArray) throw PyTrilinos::PythonException();
    PyObject * iArray = PyArray_ContiguousFromAny(indices, NPY_INT, 1, 1);
    if (!iArray)
    {
      Py_DECREF(vArray);
      throw PyTrilinos::PythonException();
    }
    int vLen = PyArray_DIM((PyArrayObject*)vArray, 0);
    int iLen = PyArray_DIM((PyArrayObject*)iArray, 0);
    if (vLen != iLen)
    {
      Py_DECREF(vArray);
      Py_DECREF(iArray);
      PyErr_Format(PyExc_ValueError,
                   "Length of values array (%d) does not equal length of "
                   "indices array (%d)",
                   vLen, iLen);
      throw PyTrilinos::PythonException();
    }
    double * vData = (double*) PyArray_DATA((PyArrayObject*)vArray);
    int *    iData = (int*)    PyArray_DATA((PyArrayObject*)iArray);
    int result = self->ReplaceMyValues(vLen, blockOffset, vData, iData);
    Py_DECREF(vArray);
    Py_DECREF(iArray);
    return result;
  }

  int SumIntoGlobalValues(PyObject * values,
                          PyObject * indices)
  {
    PyObject * vArray = PyArray_ContiguousFromAny(values, NPY_DOUBLE, 1, 1);
    if (!vArray) throw PyTrilinos::PythonException();
    PyObject * iArray = PyArray_ContiguousFromAny(indices, NPY_INT, 1, 1);
    if (!iArray)
    {
      Py_DECREF(vArray);
      throw PyTrilinos::PythonException();
    }
    int vLen = PyArray_DIM((PyArrayObject*)vArray, 0);
    int iLen = PyArray_DIM((PyArrayObject*)iArray, 0);
    if (vLen != iLen)
    {
      Py_DECREF(vArray);
      Py_DECREF(iArray);
      PyErr_Format(PyExc_ValueError,
                   "Length of values array (%d) does not equal length of "
                   "indices array (%d)",
                   vLen, iLen);
      throw PyTrilinos::PythonException();
    }
    double * vData = (double*) PyArray_DATA((PyArrayObject*)vArray);
    int *    iData = (int*)    PyArray_DATA((PyArrayObject*)iArray);
    int result = self->SumIntoGlobalValues(vLen, vData, iData);
    Py_DECREF(vArray);
    Py_DECREF(iArray);
    return result;
  }

  int SumIntoGlobalValues(int blockOffset,
                          PyObject * values,
                          PyObject * indices)
  {
    PyObject * vArray = PyArray_ContiguousFromAny(values, NPY_DOUBLE, 1, 1);
    if (!vArray) throw PyTrilinos::PythonException();
    PyObject * iArray = PyArray_ContiguousFromAny(indices, NPY_INT, 1, 1);
    if (!iArray)
    {
      Py_DECREF(vArray);
      throw PyTrilinos::PythonException();
    }
    int vLen = PyArray_DIM((PyArrayObject*)vArray, 0);
    int iLen = PyArray_DIM((PyArrayObject*)iArray, 0);
    if (vLen != iLen)
    {
      Py_DECREF(vArray);
      Py_DECREF(iArray);
      PyErr_Format(PyExc_ValueError,
                   "Length of values array (%d) does not equal length of "
                   "indices array (%d)",
                   vLen, iLen);
      throw PyTrilinos::PythonException();
    }
    double * vData = (double*) PyArray_DATA((PyArrayObject*)vArray);
    int *    iData = (int*)    PyArray_DATA((PyArrayObject*)iArray);
    int result = self->SumIntoGlobalValues(vLen, blockOffset, vData, iData);
    Py_DECREF(vArray);
    Py_DECREF(iArray);
    return result;
  }

  int SumIntoMyValues(PyObject * values,
                      PyObject * indices)
  {
    PyObject * vArray = PyArray_ContiguousFromAny(values, NPY_DOUBLE, 1, 1);
    if (!vArray) throw PyTrilinos::PythonException();
    PyObject * iArray = PyArray_ContiguousFromAny(indices, NPY_INT, 1, 1);
    if (!iArray)
    {
      Py_DECREF(vArray);
      throw PyTrilinos::PythonException();
    }
    int vLen = PyArray_DIM((PyArrayObject*)vArray, 0);
    int iLen = PyArray_DIM((PyArrayObject*)iArray, 0);
    if (vLen != iLen)
    {
      Py_DECREF(vArray);
      Py_DECREF(iArray);
      PyErr_Format(PyExc_ValueError,
                   "Length of values array (%d) does not equal length of "
                   "indices array (%d)",
                   vLen, iLen);
      throw PyTrilinos::PythonException();
    }
    double * vData = (double*) PyArray_DATA((PyArrayObject*)vArray);
    int *    iData = (int*)    PyArray_DATA((PyArrayObject*)iArray);
    int result = self->SumIntoMyValues(vLen, vData, iData);
    Py_DECREF(vArray);
    Py_DECREF(iArray);
    return result;
  }

  int SumIntoMyValues(int blockOffset,
                      PyObject * values,
                      PyObject * indices)
  {
    PyObject * vArray = PyArray_ContiguousFromAny(values, NPY_DOUBLE, 1, 1);
    if (!vArray) throw PyTrilinos::PythonException();
    PyObject * iArray = PyArray_ContiguousFromAny(indices, NPY_INT, 1, 1);
    if (!iArray)
    {
      Py_DECREF(vArray);
      throw PyTrilinos::PythonException();
    }
    int vLen = PyArray_DIM((PyArrayObject*)vArray, 0);
    int iLen = PyArray_DIM((PyArrayObject*)iArray, 0);
    if (vLen != iLen)
    {
      Py_DECREF(vArray);
      Py_DECREF(iArray);
      PyErr_Format(PyExc_ValueError,
                   "Length of values array (%d) does not equal length of "
                   "indices array (%d)",
                   vLen, iLen);
      throw PyTrilinos::PythonException();
    }
    double * vData = (double*) PyArray_DATA((PyArrayObject*)vArray);
    int *    iData = (int*)    PyArray_DATA((PyArrayObject*)iArray);
    int result = self->SumIntoMyValues(vLen, blockOffset, vData, iData);
    Py_DECREF(vArray);
    Py_DECREF(iArray);
    return result;
  }

  PyObject * __distarray__() const
  {
    return PyTrilinos::convertToDistArray(*self);
  }
}
%ignore Epetra_Vector::Epetra_Vector(Epetra_DataAccess,
                                     Epetra_BlockMap const &,
                                     double *);
%ignore Epetra_Vector::ExtractCopy(double * ) const;
%ignore Epetra_Vector::ExtractView(double **) const;
%ignore Epetra_Vector::ReplaceGlobalValues(int,double*,int*);
%ignore Epetra_Vector::ReplaceGlobalValues(int,int,double*,int*);
%ignore Epetra_Vector::ReplaceMyValues(int,double*,int*);
%ignore Epetra_Vector::ReplaceMyValues(int,int,double*,int*);
%ignore Epetra_Vector::SumIntoGlobalValues(int,double*,int*);
%ignore Epetra_Vector::SumIntoGlobalValues(int,int,double*,int*);
%ignore Epetra_Vector::SumIntoMyValues(int,double*,int*);
%ignore Epetra_Vector::SumIntoMyValues(int,int,double*,int*);
%ignore Epetra_Vector::ResetView(double *);
%include "Epetra_Vector.h"
%pythoncode
%{
  def Vector_getattr(self, name):
      if name == "array":
          a = _extractNumPyArrayFromEpetraVector(self)
          self.__dict__["array"] = a
          return a
      elif name == "shape":
          return self.array.shape
      elif name == "dtype":
          return self.array.dtype
      raise AttributeError("'%s' not an attribute of Vector" % name)
  def Vector_setattr(self, name, value):
      if name == "array":
          raise AttributeError("Cannot change Vector 'array' attribute")
      elif name == "shape":
          self.array.shape = value
      elif name == "dtype":
          raise AttributeError("Cannot change Vector 'dtype' attribute")
      else:
          self.__dict__[name] = value
  Vector.__getattr__ = Vector_getattr
  Vector.__setattr__ = Vector_setattr
  Vector.__getitem__ = lambda self, i: self.array.__getitem__(i)
  Vector.__setitem__ = lambda self, i, v: self.array.__setitem__(i,v)
  Vector.__len__     = lambda self: self.array.__len__()
  Vector.__str__     = lambda self: self.array.__str__()
  Vector.copy        = lambda self: Vector(self)
  Vector.ExtractCopy = lambda self: self.array.copy()
  Vector.ExtractView = lambda self: self.array
  class_array_add_math(Vector)
  class_array_add_comp(Vector)
%}

/////////////////////////////
// Epetra_FEVector support //
/////////////////////////////
%rename(FEVector) Epetra_FEVector;
%inline
{
  PyObject *
    _extractNumPyArrayFromEpetraFEVector(const Epetra_FEVector & source)
  {
    npy_intp dims[2] = { source.NumVectors(), source.MyLength() };
    double * data = source.Values();
    PyObject * array = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, (void*)data);
    return array;
  }
}
%extend Epetra_FEVector
{
  Epetra_FEVector(Epetra_DataAccess cv,
                  Epetra_BlockMap const & map,
                  PyObject * values,
                  int numVectors,
                  bool ignoreNonLocalEntries=false)
  {
    Epetra_FEVector * result = 0;
    int arrayDim1, arrayDim2;
    PyObject * array = PyArray_ContiguousFromAny(values, NPY_DOUBLE, 2, 2);
    if (!array) SWIG_fail;
    arrayDim1 = PyArray_DIM((PyArrayObject*)array, 0);
    arrayDim2 = PyArray_DIM((PyArrayObject*)array, 1);
    if (arrayDim1 != numVectors)
    {
      PyErr_Format(PyExc_ValueError,
                   "Array dimension %d not equal to number of vectors %d",
                   arrayDim2, numVectors);
      SWIG_fail;
    }
    result = new Epetra_FEVector(cv,
                                 map,
                                 (double*)PyArray_DATA((PyArrayObject*)array),
                                 arrayDim2,
                                 arrayDim1,
                                 ignoreNonLocalEntries);
    Py_DECREF(array);
    return result;
  fail:
    Py_XDECREF(array);
    return NULL;
  }

  int ReplaceGlobalValues(PyObject * gids,
                          PyObject * values,
                          int vectorId=0)
  {
    PyObject * gid_array = 0;
    PyObject * val_array = 0;
    int numGid, numVal;
    int err = 0;
    gid_array = PyArray_ContiguousFromAny(gids  , NPY_INT   , 1, 1);
    if (!gid_array) SWIG_fail;
    val_array = PyArray_ContiguousFromAny(values, NPY_DOUBLE, 1, 1);
    if (!val_array) SWIG_fail;
    numGid = PyArray_DIM((PyArrayObject*)gid_array, 0);
    numVal = PyArray_DIM((PyArrayObject*)val_array, 0);
    if (numGid != numVal)
    {
      PyErr_Format(PyExc_ValueError,
                   "Number of GIDs %d != number of values %d",
                   numGid, numVal);
      SWIG_fail;
    }
    err =
      self->ReplaceGlobalValues(numGid,
                                (int*)PyArray_DATA((PyArrayObject*)gid_array),
                                (double*)PyArray_DATA((PyArrayObject*)val_array),
                                vectorId);
    Py_DECREF(gid_array);
    Py_DECREF(val_array);
    return err;
  fail:
    Py_XDECREF(gid_array);
    Py_XDECREF(val_array);
    return -1;
  }

  int SumIntoGlobalValues(PyObject * gids,
                          PyObject * values,
                          int vectorId=0)
  {
    PyObject * gid_array = 0;
    PyObject * val_array = 0;
    int numGid, numVal;
    int err = 0;
    gid_array = PyArray_ContiguousFromAny(gids  , NPY_INT   , 1, 1);
    if (!gid_array) SWIG_fail;
    val_array = PyArray_ContiguousFromAny(values, NPY_DOUBLE, 1, 1);
    if (!val_array) SWIG_fail;
    numGid = PyArray_DIM((PyArrayObject*)gid_array, 0);
    numVal = PyArray_DIM((PyArrayObject*)val_array, 0);
    if (numGid != numVal)
    {
      PyErr_Format(PyExc_ValueError,
                   "Number of GIDs %d != number of values %d",
                   numGid, numVal);
      SWIG_fail;
    }
    err =
      self->SumIntoGlobalValues(numGid,
                                (int*)PyArray_DATA((PyArrayObject*)gid_array),
                                (double*)PyArray_DATA((PyArrayObject*)val_array),
                                vectorId);
    Py_DECREF(gid_array);
    Py_DECREF(val_array);
    return err;
  fail:
    Py_XDECREF(gid_array);
    Py_XDECREF(val_array);
    return -1;
  }

  PyObject * __distarray__() const
  {
    return PyTrilinos::convertToDistArray(*self);
  }
}
%ignore Epetra_FEVector::Epetra_FEVector(Epetra_DataAccess,
                                         const Epetra_BlockMap &,
                                         double *,
                                         int,
                                         int,
                                         bool);
%ignore Epetra_FEVector::Epetra_FEVector(Epetra_DataAccess,
                                         const Epetra_BlockMap &,
                                         double **,
                                         int,
                                         int,
                                         bool);
%ignore Epetra_FEVector::ReplaceGlobalValues(int,int*,double*);
%ignore Epetra_FEVector::SumIntoGlobalValues(int,int*,double*);
%include "Epetra_FEVector.h"
%pythoncode
%{
  def FEVector_getattr(self, name):
      if name == "array":
          a = _extractNumPyArrayFromEpetraFEVector(self)
          self.__dict__["array"] = a
          return a
      elif name == "shape":
          return self.array.shape
      elif name == "dtype":
          return self.array.dtype
      raise AttributeError("'%s' not an attribute of FEVector" % name)
  def FEVector_setattr(self, name, value):
      if name == "array":
          raise AttributeError("Cannot change FEVector 'array' attribute")
      elif name == "shape":
          self.array.shape = value
      elif name == "dtype":
          raise AttributeError("Cannot change FEVector 'dtype' attribute")
      else:
          self.__dict__[name] = value
  FEVector.__getattr__ = FEVector_getattr
  FEVector.__setattr__ = FEVector_setattr
  FEVector.__getitem__ = lambda self, i: self.array.__getitem__(i)
  FEVector.__setitem__ = lambda self, i, v: self.array.__setitem__(i, v)
  FEVector.__len__     = lambda self: self.array.__len__()
  FEVector.__str__     = lambda self: self.array.__str__()
  FEVector.copy        = lambda self: Vector(self)
  FEVector.ExtractCopy = lambda self: self.array.copy()
  FEVector.ExtractView = lambda self: self.array
  class_array_add_math(FEVector)
  class_array_add_comp(FEVector)
%}

// Access RCP methods
%inline
{
  int strong_count(Teuchos::RCP< Epetra_Vector > & ev)
  {
    return ev.strong_count();
  }
  int weak_count(Teuchos::RCP< Epetra_Vector > & ev)
  {
    return ev.weak_count();
  }
  int total_count(Teuchos::RCP< Epetra_Vector > & ev)
  {
    return ev.total_count();
  }
}
