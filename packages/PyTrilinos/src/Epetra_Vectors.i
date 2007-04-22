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
// Epetra includes
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntVector.h"

// Local includes
#include "Epetra_NumPyIntVector.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyFEVector.h"
%}

//////////////
// Typemaps //
//////////////
%epetra_array_output_typemaps(Epetra_IntVector,   Epetra_NumPyIntVector  )
%epetra_array_output_typemaps(Epetra_MultiVector, Epetra_NumPyMultiVector)
%epetra_array_output_typemaps(Epetra_Vector,      Epetra_NumPyVector     )
%epetra_array_output_typemaps(Epetra_FEVector,    Epetra_NumPyFEVector   )
%epetra_array_argout_typemaps(IntVector)
%epetra_array_argout_typemaps(MultiVector)
%epetra_array_argout_typemaps(Vector)
%epetra_array_argout_typemaps(FEVector)
%epetra_array_director_typemaps(MultiVector)
%epetra_array_director_typemaps(Vector)

//////////////////////////////
// Epetra_IntVector support //
//////////////////////////////
%inline {struct IntVector{ };}
%include "Epetra_IntVector.h"

////////////////////////////////
// Epetra_MultiVector support //
////////////////////////////////
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
%inline {struct MultiVector{ };}
%include "Epetra_MultiVector.h"

///////////////////////////
// Epetra_Vector support //
///////////////////////////
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
%inline {struct Vector{ };}
%include "Epetra_Vector.h"

/////////////////////////////
// Epetra_FEVector support //
/////////////////////////////
%ignore Epetra_FEVector::ReplaceGlobalValues(int,int*,double*);
%ignore Epetra_FEVector::SumIntoGlobalValues(int,int*,double*);
%inline {struct FEVector{ };}
%include "Epetra_FEVector.h"

///////////////////////////////////
// Epetra_NumPyIntVector support //
///////////////////////////////////
%rename(NumPyIntVector) Epetra_NumPyIntVector;
%epetra_numpy_ctor_exception(Epetra_NumPyIntVector)
%include "Epetra_NumPyIntVector.h"
%pythoncode %{
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
        return numpy.less(self.array,other)
    def __le__(self,other):
        return numpy.less_equal(self.array,other)
    def __eq__(self,other):
        return numpy.equal(self.array,other)
    def __ne__(self,other):
        return numpy.not_equal(self.array,other)
    def __gt__(self,other):
        return numpy.greater(self.array,other)
    def __ge__(self,other):
        return numpy.greater_equal(self.array,other)
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

/////////////////////////////////////
// Epetra_NumPyMultiVector support //
/////////////////////////////////////
%rename(NumPyMultiVector) Epetra_NumPyMultiVector;
%epetra_numpy_ctor_exception(Epetra_NumPyMultiVector)
%include "Epetra_NumPyMultiVector.h"
%pythoncode %{
class MultiVector(UserArray,NumPyMultiVector):
    def __init__(self, *args):
        """
        __init__(self, BlockMap map, int numVectors, bool zeroOut=True) -> MultiVector
        __init__(self, MultiVector source) -> MultiVector
        __init__(self, BlockMap map, PyObject array) -> MultiVector
        __init__(self, DataAccess CV, MultiVector source) -> MultiVector
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
        return numpy.less(self.array,other)
    def __le__(self,other):
        return numpy.less_equal(self.array,other)
    def __eq__(self,other):
        return numpy.equal(self.array,other)
    def __ne__(self,other):
        return numpy.not_equal(self.array,other)
    def __gt__(self,other):
        return numpy.greater(self.array,other)
    def __ge__(self,other):
        return numpy.greater_equal(self.array,other)
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
%}

////////////////////////////////
// Epetra_NumPyVector support //
////////////////////////////////
%rename(NumPyVector) Epetra_NumPyVector;
%epetra_numpy_ctor_exception(Epetra_NumPyVector)
%include "Epetra_NumPyVector.h"
%pythoncode %{
class Vector(UserArray,NumPyVector):
    def __init__(self, *args):
        """
        __init__(self, BlockMap map, bool zeroOut=True) -> Vector
        __init__(self, Vector source) -> Vector
        __init__(self, BlockMap map, PyObject array) -> Vector
        __init__(self, DataAccess CV, Vector source) -> Vector
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
        return numpy.less(self.array,other)
    def __le__(self,other):
        return numpy.less_equal(self.array,other)
    def __eq__(self,other):
        return numpy.equal(self.array,other)
    def __ne__(self,other):
        return numpy.not_equal(self.array,other)
    def __gt__(self,other):
        return numpy.greater(self.array,other)
    def __ge__(self,other):
        return numpy.greater_equal(self.array,other)
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
%}

////////////////////////////////
// Epetra_NumPyFEVector support //
////////////////////////////////
%rename(NumPyFEVector) Epetra_NumPyFEVector;
%epetra_numpy_ctor_exception(Epetra_NumPyFEVector)
%include "Epetra_NumPyFEVector.h"
%pythoncode %{
class FEVector(UserArray,NumPyFEVector):
    def __init__(self, *args):
        """
        __init__(self, BlockMap map, bool zeroOut=True) -> FEVector
        __init__(self, FEVector source) -> FEVector
        __init__(self, BlockMap map, PyObject array) -> FEVector
        __init__(self, DataAccess CV, Vector source) -> FEVector
        __init__(self, DataAccess CV, MultiVector source, PyObject index) -> FEVector
        __init__(self, PyObject array) -> FEVector
        """
        NumPyFEVector.__init__(self, *args)
        self.__initArray__()
    def __initArray__(self):
        UserArray.__init__(self, self.ExtractView(), dtype="d", copy=False)
    def __str__(self):
        return str(self.array)
    def __lt__(self,other):
        return numpy.less(self.array,other)
    def __le__(self,other):
        return numpy.less_equal(self.array,other)
    def __eq__(self,other):
        return numpy.equal(self.array,other)
    def __ne__(self,other):
        return numpy.not_equal(self.array,other)
    def __gt__(self,other):
        return numpy.greater(self.array,other)
    def __ge__(self,other):
        return numpy.greater_equal(self.array,other)
    def __getattr__(self, key):
        # This should get called when the FEVector is accessed after not properly
        # being initialized
        if not "array" in self.__dict__:
            self.__initArray__()
        try:
            return self.array.__getattribute__(key)
        except AttributeError:
            return FEVector.__getattribute__(self, key)
    def __setattr__(self, key, value):
        "Handle 'this' properly and protect the 'array' attribute"
        if key == "this":
            NumPyFEVector.__setattr__(self, key, value)
        else:
            if key == "array":
                if key in self.__dict__:
                    raise AttributeError, "Cannot change Epetra.FEVector array attribute"
            UserArray.__setattr__(self, key, value)
_Epetra.NumPyFEVector_swigregister(FEVector)
%}
