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
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseOperator.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialSymDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialDenseSVD.h"

// Local interface includes
#include "Epetra_NumPySerialDenseMatrix.h"
#include "Epetra_NumPySerialDenseVector.h"
#include "Epetra_NumPyIntSerialDenseVector.h"
%}

// Ignore directives
%ignore Epetra_SerialDenseMatrix::operator()(int,int) const;
%ignore Epetra_SerialDenseMatrix::A() const;
%ignore Epetra_SerialDenseVector::operator()(int);
%ignore Epetra_SerialDenseVector::operator()(int) const;
%ignore Epetra_IntSerialDenseMatrix::operator()(int,int) const;
%ignore Epetra_IntSerialDenseMatrix::A() const;
%ignore Epetra_IntSerialDenseVector::operator()(int);
%ignore Epetra_IntSerialDenseVector::operator()(int) const;

// Rename directives
%rename(IntSerialDenseMatrix     ) Epetra_IntSerialDenseMatrix;
%rename(NumPyIntSerialDenseVector) Epetra_NumPyIntSerialDenseVector;
%rename(SerialDenseOperator   	 ) Epetra_SerialDenseOperator;
%rename(NumPySerialDenseMatrix	 ) Epetra_NumPySerialDenseMatrix;
%rename(SerialSymDenseMatrix  	 ) Epetra_SerialSymDenseMatrix;
%rename(NumPySerialDenseVector	 ) Epetra_NumPySerialDenseVector;
%rename(SerialDenseSolver     	 ) Epetra_SerialDenseSolver;
%rename(SerialDenseSVD        	 ) Epetra_SerialDenseSVD;

// These are place-holders on the C++ side for the python-defined
// SerialDenseMatrix, SerialDenseVector and IntSerialDenseVector
%inline {
  struct SerialDenseMatrix   { };
  struct SerialDenseVector   { };
  struct IntSerialDenseVector{ };
}

// Typemaps
%define TYPEMAP_OUT(vector,numPyVector)
%typemap(out) vector * {
  numPyVector * vec = new numPyVector(*$1);
  static swig_type_info *ty = SWIG_TypeQuery("numPyVector *");
  $result = SWIG_NewPointerObj(vec, ty, 1);
}
%enddef

TYPEMAP_OUT(Epetra_SerialDenseMatrix,   Epetra_NumPySerialDenseMatrix   )
TYPEMAP_OUT(Epetra_SerialDenseVector,   Epetra_NumPySerialDenseVector   )
TYPEMAP_OUT(Epetra_IntSerialDenseVector,Epetra_NumPyIntSerialDenseVector)

// Epetra include directives
%include "Epetra_IntSerialDenseMatrix.h"
%include "Epetra_IntSerialDenseVector.h"
%include "Epetra_SerialDenseOperator.h"
%include "Epetra_SerialDenseMatrix.h"
%include "Epetra_SerialSymDenseMatrix.h"
%include "Epetra_SerialDenseVector.h"
%include "Epetra_SerialDenseSolver.h"
%include "Epetra_SerialDenseSVD.h"

// Local interface include directives
%include "Epetra_NumPySerialDenseMatrix.h"
%include "Epetra_NumPySerialDenseVector.h"
%include "Epetra_NumPyIntSerialDenseVector.h"

// Python code
%pythoncode %{

from UserArray import *

class SerialDenseMatrix(UserArray,NumPySerialDenseMatrix):
    def __init__(self, *args):
      	"""
      	__init__(self, bool set_object_label=True) -> SerialDenseMatrix
      	__init__(self, int numRows, int numCols, bool set_object_label=True) -> SerialDenseMatrix
      	__init__(self, PyObject array, bool set_object_label=True) -> SerialDenseMatrix
      	__init__(self, SerialDenseMatrix source) -> SerialDenseMatrix
      	"""
        NumPySerialDenseMatrix.__init__(self, *args)
        self.__initArray__()
    def __initArray__(self):
        UserArray.__init__(self,self.A(),'d',copy=False,savespace=True)
        self.__protected = True
    def __str__(self):
        return str(self.array)
    def __getattr__(self, key):
        # This should get called when the SerialDenseMatrix is accessed after
        # not properly being initialized
        if not 'array' in self.__dict__:
            self.__initArray__()
        return self.__dict__[key]
    def __setattr__(self, key, value):
        "Protect the 'array' and 'shape' attributes"
        if key in self.__dict__:
            if self.__protected:
                if key == "array":
                    raise AttributeError, "Cannot change Epetra.SerialDenseMatrix array attribute"
                if key == "shape":
                    raise AttributeError, "Cannot change Epetra.SerialDenseMatrix shape attribute"
        UserArray.__setattr__(self, key, value)
    def __call__(self,i,j):
        "__call__(self, int i, int j) -> double"
        return self.__getitem__(i,j)
    def Shape(self,numRows,numCols):
        "Shape(self, int numRows, int numCols) -> int"
        result = NumPySerialDenseMatrix.Shape(self,numRows,numCols)
        self.__protected = False
        UserArray.__init__(self,self.A(),'d',copy=False,savespace=True)
        self.__protected = True
    def Reshape(self,numRows,numCols):
        "Reshape(self, int numRows, int numCols) -> int"
        result = NumPySerialDenseMatrix.Reshape(self,numRows,numCols)
        self.__protected = False
        UserArray.__init__(self,self.A(),'d',copy=False,savespace=True)
        self.__protected = True
_Epetra.NumPySerialDenseMatrix_swigregister(SerialDenseMatrix)

class SerialDenseVector(UserArray,NumPySerialDenseVector):
    def __init__(self, *args):
      	"""
      	__init__(self) -> SerialDenseVector
      	__init__(self, int length) -> SerialDenseVector
      	__init__(self, PyObject array) -> SerialDenseVector
      	__init__(self, SerialDenseVector source) -> SerialDenseVector
      	"""
        NumPySerialDenseVector.__init__(self, *args)
        self.CheckForError()
        self.__initArray__()
    def __initArray__(self):
        UserArray.__init__(self,self.Values(),'d',copy=False,savespace=True)
        self.__protected = True
    def __str__(self):
        return str(self.array)
    def __getattr__(self, key):
        # This should get called when the SerialDenseVector is accessed after
        # not properly being initialized
        if not 'array' in self.__dict__:
            self.__initArray__()
        return self.__dict__[key]
    def __setattr__(self, key, value):
        "Protect the 'array' attribute"
        if key in self.__dict__:
            if self.__protected:
                if key == "array":
                    raise AttributeError, "Cannot change Epetra.SerialDenseVector array attribute"
        UserArray.__setattr__(self, key, value)
    def __call__(self,i):
        "__call__(self, int i) -> double"
        return self.__getitem__(i)
    def Size(self,length):
        "Size(self, int length) -> int"
        result = NumPySerialDenseVector.Size(self,length)
        self.__protected = False
        UserArray.__init__(self,self.Values(),'d',copy=False,savespace=True)
        self.__protected = True
    def Resize(self,length):
        "Resize(self, int length) -> int"
        result = NumPySerialDenseVector.Resize(self,length)
        self.__protected = False
        UserArray.__init__(self,self.Values(),'d',copy=False,savespace=True)
        self.__protected = True
_Epetra.NumPySerialDenseVector_swigregister(SerialDenseVector)

class IntSerialDenseVector(UserArray,NumPyIntSerialDenseVector):
    def __init__(self, *args):
      	"""
      	__init__(self) -> IntSerialDenseVector
      	__init__(self, int length) -> IntSerialDenseVector
      	__init__(self, PyObject array) -> IntSerialDenseVector
      	__init__(self, IntSerialDenseVector source) -> IntSerialDenseVector
      	"""
        NumPyIntSerialDenseVector.__init__(self, *args)
        self.CheckForError()
        self.__initArray__()
    def __initArray__(self):
        UserArray.__init__(self,self.Values(),'i',copy=False,savespace=True)
        self.__protected = True
    def __str__(self):
        return str(self.array)
    def __getattr__(self, key):
        # This should get called when the IntSerialDenseVector is accessed after
        # not properly being initialized
        if not 'array' in self.__dict__:
            self.__initArray__()
        return self.__dict__[key]
    def __setattr__(self, key, value):
        "Protect the 'array' attribute"
        if key in self.__dict__:
            if self.__protected:
                if key == "array":
                    raise AttributeError, "Cannot change Epetra.IntSerialDenseVector array attribute"
        UserArray.__setattr__(self, key, value)
    def __call__(self,i):
        "__call__(self, int i) -> int"
        return self.__getitem__(i)
    def Size(self,length):
        "Size(self, int length) -> int"
        result = NumPyIntSerialDenseVector.Size(self,length)
        self.__protected = False
        UserArray.__init__(self,self.Values(),'i',copy=False,savespace=True)
        self.__protected = True
    def Resize(self,length):
        "Resize(self, int length) -> int"
        result = NumPyIntSerialDenseVector.Resize(self,length)
        self.__protected = False
        UserArray.__init__(self,self.Values(),'i',copy=False,savespace=True)
        self.__protected = True
_Epetra.NumPyIntSerialDenseVector_swigregister(IntSerialDenseVector)

%}

// Epetra_SerialSpdDenseSolver is apparently not built
//#include "Epetra_SerialSpdDenseSolver.h"
//%rename(SerialSpdDenseSolver  ) Epetra_SerialSpdDenseSolver;
//%include "Epetra_SerialSpdDenseSolver.h"
