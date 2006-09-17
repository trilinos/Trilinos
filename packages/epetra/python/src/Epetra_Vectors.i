// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//            PyTrilinos.Epetra: Python Interface to Epetra
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

%{
// Epetra includes
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntVector.h"

// Local includes
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyIntVector.h"
%}

// Ignore directives
%ignore Epetra_MultiVector::operator()(int) const;
// The following methods are given functionality in the derived class
// Epetra_NumPyMultiVector
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
// The following are expert methods not supported in python
%ignore Epetra_MultiVector::ResetView(double **);
%ignore Epetra_MultiVector::Pointers() const;
// The following methods are given functionality in the derived class
// Epetra_NumPyVector
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
// The following is an expert method not supported in python
%ignore Epetra_Vector::ResetView(double *);

// Rename directives
%rename(NumPyMultiVector) Epetra_NumPyMultiVector;
%rename(NumPyVector     ) Epetra_NumPyVector;
%rename(NumPyIntVector  ) Epetra_NumPyIntVector;
%rename(FEVector        ) Epetra_FEVector;

// These are place-holders on the C++ side for the python-defined
// MultiVector, Vector and IntVector classes
%inline {
  struct MultiVector { };
  struct Vector      { };
  struct IntVector   { };
}

// Exceptions.  The Epetra_NumPy*Vector class constructors can raise
// python exceptions, but swig must be told explicitely to look for
// them.
NUMPY_CONSTRUCTOR_EXCEPTION_HANDLER(Epetra_NumPyMultiVector)
NUMPY_CONSTRUCTOR_EXCEPTION_HANDLER(Epetra_NumPyVector     )
NUMPY_CONSTRUCTOR_EXCEPTION_HANDLER(Epetra_NumPyIntVector  )

// Typemaps
TYPEMAP_OUT(Epetra_MultiVector,Epetra_NumPyMultiVector)
TYPEMAP_OUT(Epetra_Vector,     Epetra_NumPyVector     )
TYPEMAP_OUT(Epetra_IntVector,  Epetra_NumPyIntVector  )

// Include directives for Epetra
%include "Epetra_IntVector.h"
%include "Epetra_MultiVector.h"
%include "Epetra_Vector.h"
%include "Epetra_FEVector.h"

// Local interface includes
%include "Epetra_NumPyMultiVector.h"
%include "Epetra_NumPyVector.h"
%include "Epetra_NumPyIntVector.h"

// Python code.  Here we define the Epetra.MultiVector, Epetra.Vector
// and Epetra.IntVector python classes, which multiply-inherit from
// the numpy UserArray class (making these classes numpy arrays) and
// the Epetra_NumPyMultiVector, Epetra_NumPyVector or
// Epetra_NumPyIntVector class (making these classes also Epetra
// objects).
%pythoncode %{

class MultiVector(UserArray,NumPyMultiVector):
    def __init__(self, *args):
        """
        __init__(self, BlockMap map, int numVectors, bool zeroOut=True) -> MultiVector
        __init__(self, MultiVector source) -> MultiVector
        __init__(self, BlockMap map, PyObject array) -> MultiVector
        __init__(self, DataAccess CV, MultiVector source, PyObject range) -> MultiVector
        __init__(self, PyObject array) -> MultiVector
        """
        NumPyMultiVector.__init__(self, *args)
        self.__initArray__()
    def __initArray__(self):
        UserArray.__init__(self, self.ExtractView(), dtype="d", copy=False)
    def __str__(self):
        return str(self.array)
    def __lt__(self,other):
        return less(self.array,other)
    def __le__(self,other):
        return less_equal(self.array,other)
    def __eq__(self,other):
        return equal(self.array,other)
    def __ne__(self,other):
        return not_equal(self.array,other)
    def __gt__(self,other):
        return greater(self.array,other)
    def __ge__(self,other):
        return greater_equal(self.array,other)
    def __getattr__(self, key):
        # This should get called when the MultiVector is accessed after not
        # properly being initialized
        if not "array" in self.__dict__:
            self.__initArray__()
        try:
            return self.array.__getattribute__(key)
        except AttributeError:
            return MultiVector.__getattribute__(self, key)
    def __setattr__(self, key, value):
        "Handle 'this' properly and protect the 'array' and 'shape' attributes"
        if key == "this":
            NumPyMultiVector.__setattr__(self, key, value)
        else:
            if key == "array":
                if key in self.__dict__:
                    raise AttributeError, "Cannot change Epetra.MultiVector array attribute"
            elif key == "shape":
                value = tuple(value)
                if len(value) < 2:
                    raise ValueError, "Epetra.MultiVector shape is " + str(value) + \
                      " but must have minimum of 2 elements"
            UserArray.__setattr__(self, key, value)
_Epetra.NumPyMultiVector_swigregister(MultiVector)

class Vector(UserArray,NumPyVector):
    def __init__(self, *args):
        """
        __init__(self, BlockMap map, bool zeroOut=True) -> Vector
        __init__(self, Vector source) -> Vector
        __init__(self, BlockMap map, PyObject array) -> Vector
        __init__(self, DataAccess CV, MultiVector source, PyObject index) -> Vector
        __init__(self, PyObject array) -> Vector
        """
        NumPyVector.__init__(self, *args)
        self.__initArray__()
    def __initArray__(self):
        UserArray.__init__(self, self.ExtractView(), dtype="d", copy=False)
    def __str__(self):
        return str(self.array)
    def __lt__(self,other):
        return less(self.array,other)
    def __le__(self,other):
        return less_equal(self.array,other)
    def __eq__(self,other):
        return equal(self.array,other)
    def __ne__(self,other):
        return not_equal(self.array,other)
    def __gt__(self,other):
        return greater(self.array,other)
    def __ge__(self,other):
        return greater_equal(self.array,other)
    def __getattr__(self, key):
        # This should get called when the Vector is accessed after not properly
        # being initialized
        if not "array" in self.__dict__:
            self.__initArray__()
        try:
            return self.array.__getattribute__(key)
        except AttributeError:
            return Vector.__getattribute__(self, key)
    def __setattr__(self, key, value):
        "Handle 'this' properly and protect the 'array' attribute"
        if key == "this":
            NumPyVector.__setattr__(self, key, value)
        else:
            if key == "array":
                if key in self.__dict__:
                    raise AttributeError, "Cannot change Epetra.Vector array attribute"
            UserArray.__setattr__(self, key, value)
_Epetra.NumPyVector_swigregister(Vector)

class IntVector(UserArray,NumPyIntVector):
    def __init__(self, *args):
        """
        __init__(self, BlockMap map, bool zeroOut=True) -> IntVector
        __init__(self, IntVector source) -> IntVector
        __init__(self, BlockMap map, PyObject array) -> IntVector
        __init__(self, PyObject array) -> IntVector
        """
        NumPyIntVector.__init__(self, *args)
        self.__initArray__()
    def __initArray__(self):
        UserArray.__init__(self, self.ExtractView(), dtype="i", copy=False)
    def __str__(self):
        return str(self.array)
    def __lt__(self,other):
        return less(self.array,other)
    def __le__(self,other):
        return less_equal(self.array,other)
    def __eq__(self,other):
        return equal(self.array,other)
    def __ne__(self,other):
        return not_equal(self.array,other)
    def __gt__(self,other):
        return greater(self.array,other)
    def __ge__(self,other):
        return greater_equal(self.array,other)
    def __getattr__(self, key):
        # This should get called when the IntVector is accessed after not
        # properly being initialized
        if not "array" in self.__dict__:
            self.__initArray__()
        try:
            return self.array.__getattribute__(key)
        except AttributeError:
            return IntVector.__getattribute__(self, key)
    def __setattr__(self, key, value):
        "Handle 'this' properly and protect the 'array' attribute"
        if key == "this":
            NumPyIntVector.__setattr__(self, key, value)
        else:
            if key == "array":
                if key in self.__dict__:
                    raise AttributeError, "Cannot change Epetra.IntVector array attribute"
            UserArray.__setattr__(self, key, value)
_Epetra.NumPyIntVector_swigregister(IntVector)

%}
