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
// Local interface includes
#include "Epetra_NumPyIntSerialDenseMatrix.h"
#include "Epetra_NumPyIntSerialDenseVector.h"
#include "Epetra_NumPySerialDenseMatrix.h"
#include "Epetra_NumPySerialSymDenseMatrix.h"
#include "Epetra_NumPySerialDenseVector.h"

// Epetra includes
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
// wrappers.  However, here it appears to causes errors, so we turn it
// off.
%feature("compactdefaultargs", "0");

////////////////////////////////////////////////
// Typemaps for Teuchos::RCP< Epetra_NumPy* > //
////////////////////////////////////////////////
%teuchos_rcp_epetra_numpy(IntSerialDenseMatrix)
%teuchos_rcp_epetra_numpy(IntSerialDenseVector)
%teuchos_rcp_epetra_numpy(SerialDenseMatrix   )
%teuchos_rcp_epetra_numpy(SerialSymDenseMatrix)
%teuchos_rcp_epetra_numpy(SerialDenseVector   )

/////////////////////////////////////////
// Epetra_IntSerialDenseMatrix support //
/////////////////////////////////////////
%ignore Epetra_IntSerialDenseMatrix::operator()(int,int) const;
%ignore Epetra_IntSerialDenseMatrix::A() const;
%ignore Epetra_IntSerialDenseMatrix::MakeViewOf;
%inline
{
  struct IntSerialDenseMatrix{ };
}
%include "Epetra_IntSerialDenseMatrix.h"

//////////////////////////////////////////////
// Epetra_NumPyIntSerialDenseMatrix support //
//////////////////////////////////////////////
%rename(NumPyIntSerialDenseMatrix) Epetra_NumPyIntSerialDenseMatrix;
%include "Epetra_NumPyIntSerialDenseMatrix.h"
%pythoncode
%{
class IntSerialDenseMatrix(UserArray,NumPyIntSerialDenseMatrix):
    def __init__(self, *args):
      	"""
      	__init__(self) -> IntSerialDenseMatrix
      	__init__(self, int numRows, int numCols) -> IntSerialDenseMatrix
      	__init__(self, PyObject array) -> IntSerialDenseMatrix
      	__init__(self, IntSerialDenseMatrix source) -> IntSerialDenseMatrix
      	"""
        NumPyIntSerialDenseMatrix.__init__(self, *args)
        self.__initArray__()
    def __initArray__(self):
        self.array = self.A()
        self.__protected = True
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
        # This should get called when the IntSerialDenseMatrix is accessed after
        # not properly being initialized
        if not "array" in self.__dict__:
            self.__initArray__()
        try:
            return self.array.__getattribute__(key)
        except AttributeError:
            return IntSerialDenseMatrix.__getattribute__(self, key)
    def __setattr__(self, key, value):
        "Handle 'this' attribute properly and protect the 'array' and 'shape' attributes"
        if key == "this":
            NumPyIntSerialDenseMatrix.__setattr__(self, key, value)
        else:
            if key in self.__dict__:
                if self.__protected:
                    if key == "array":
                        raise AttributeError, \
                              "Cannot change Epetra.IntSerialDenseMatrix array attribute"
                    if key == "shape":
                        raise AttributeError, \
                              "Cannot change Epetra.IntSerialDenseMatrix shape attribute"
            UserArray.__setattr__(self, key, value)
    def __getitem__(self, index):
        """
        __getitem__(self,int,int) -> int
        __getitem__(self,int,slice) -> array
        __getitem__(self,slice,int) -> array
        __getitem__(self,slice,slice) -> array
        """
        return self.array[index]
    def Shape(self,numRows,numCols):
        "Shape(self, int numRows, int numCols) -> int"
        result = NumPyIntSerialDenseMatrix.Shape(self,numRows,numCols)
        self.__protected = False
        self.__initArray__()
        return result
    def Reshape(self,numRows,numCols):
        "Reshape(self, int numRows, int numCols) -> int"
        result = NumPyIntSerialDenseMatrix.Reshape(self,numRows,numCols)
        self.__protected = False
        self.__initArray__()
        return result
_Epetra.NumPyIntSerialDenseMatrix_swigregister(IntSerialDenseMatrix)
%}

/////////////////////////////////////////
// Epetra_IntSerialDenseVector support //
/////////////////////////////////////////
%ignore Epetra_IntSerialDenseVector::operator()(int);
%ignore Epetra_IntSerialDenseVector::operator()(int) const;
%inline
{
  struct IntSerialDenseVector{ };
}
%include "Epetra_IntSerialDenseVector.h"

//////////////////////////////////////////////
// Epetra_NumPyIntSerialDenseVector support //
//////////////////////////////////////////////
%rename(NumPyIntSerialDenseVector) Epetra_NumPyIntSerialDenseVector;
%include "Epetra_NumPyIntSerialDenseVector.h"
%pythoncode
%{
class IntSerialDenseVector(UserArray,NumPyIntSerialDenseVector):
    def __init__(self, *args):
      	"""
      	__init__(self) -> IntSerialDenseVector
      	__init__(self, int length) -> IntSerialDenseVector
      	__init__(self, PyObject array) -> IntSerialDenseVector
      	__init__(self, IntSerialDenseVector source) -> IntSerialDenseVector
      	"""
        NumPyIntSerialDenseVector.__init__(self, *args)
        self.__initArray__()
    def __initArray__(self):
        self.array = self.Values()
        self.__protected = True
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
        # This should get called when the IntSerialDenseVector is accessed after
        # not properly being initialized
        if not "array" in self.__dict__:
            self.__initArray__()
        try:
            return self.array.__getattribute__(key)
        except AttributeError:
            return IntSerialDenseVector.__getattribute__(self, key)
    def __setattr__(self, key, value):
        "Handle 'this' attribute properly and protect the 'array' attribute"
        if key == "this":
            NumPyIntSerialDenseVector.__setattr__(self, key, value)
        else:
            if key in self.__dict__:
                if self.__protected:
                    if key == "array":
                        raise AttributeError, \
                              "Cannot change Epetra.IntSerialDenseVector array attribute"
            UserArray.__setattr__(self, key, value)
    def __call__(self,i):
        "__call__(self, int i) -> int"
        return self.__getitem__(i)
    def Size(self,length):
        "Size(self, int length) -> int"
        result = NumPyIntSerialDenseVector.Size(self,length)
        self.__protected = False
        self.__initArray__()
        return result
    def Resize(self,length):
        "Resize(self, int length) -> int"
        result = NumPyIntSerialDenseVector.Resize(self,length)
        self.__protected = False
        self.__initArray__()
        return result
_Epetra.NumPyIntSerialDenseVector_swigregister(IntSerialDenseVector)
%}

////////////////////////////////////////
// Epetra_SerialDenseOperator support //
////////////////////////////////////////
%rename(SerialDenseOperator) Epetra_SerialDenseOperator;
%include "Epetra_SerialDenseOperator.h"

//////////////////////////////////////
// Epetra_SerialDenseMatrix support //
//////////////////////////////////////
%ignore Epetra_SerialDenseMatrix::operator()(int,int) const;
%ignore Epetra_SerialDenseMatrix::A() const;
%inline
{
  struct SerialDenseMatrix{ };
}
%include "Epetra_SerialDenseMatrix.h"

///////////////////////////////////////////
// Epetra_NumPySerialDenseMatrix support //
///////////////////////////////////////////
%rename(NumPySerialDenseMatrix) Epetra_NumPySerialDenseMatrix;
%include "Epetra_NumPySerialDenseMatrix.h"
%pythoncode
%{
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
        self.array = self.A()
        self.__protected = True
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
        # This should get called when the SerialDenseMatrix is accessed after
        # not properly being initialized
        if not "array" in self.__dict__:
            self.__initArray__()
        try:
            return self.array.__getattribute__(key)
        except AttributeError:
            return SerialDenseMatrix.__getattribute__(self, key)
    def __setattr__(self, key, value):
        "Handle 'this' attribute properly and protect the 'array' and 'shape' attributes"
        if key == "this":
            NumPySerialDenseMatrix.__setattr__(self, key, value)
        else:
            if key in self.__dict__:
                if self.__protected:
                    if key == "array":
                        raise AttributeError, \
                              "Cannot change Epetra.SerialDenseMatrix array attribute"
                    if key == "shape":
                        raise AttributeError, \
                              "Cannot change Epetra.SerialDenseMatrix shape attribute"
            UserArray.__setattr__(self, key, value)
    def __getitem__(self, index):
        """
        __getitem__(self,int,int) -> int
        __getitem__(self,int,slice) -> array
        __getitem__(self,slice,int) -> array
        __getitem__(self,slice,slice) -> array
        """
        return self.array[index]
    def Shape(self,numRows,numCols):
        "Shape(self, int numRows, int numCols) -> int"
        result = NumPySerialDenseMatrix.Shape(self,numRows,numCols)
        self.__protected = False
        self.__initArray__()
        return result
    def Reshape(self,numRows,numCols):
        "Reshape(self, int numRows, int numCols) -> int"
        result = NumPySerialDenseMatrix.Reshape(self,numRows,numCols)
        self.__protected = False
        self.__initArray__()
        return result
_Epetra.NumPySerialDenseMatrix_swigregister(SerialDenseMatrix)
%}

/////////////////////////////////////////
// Epetra_SerialSymDenseMatrix support //
/////////////////////////////////////////
%inline
{
  struct SerialSymDenseMatrix{ };
}
%include "Epetra_SerialSymDenseMatrix.h"

///////////////////////////////////////////
// Epetra_NumPySerialSymDenseMatrix support //
///////////////////////////////////////////
%rename(NumPySerialSymDenseMatrix) Epetra_NumPySerialSymDenseMatrix;
%include "Epetra_NumPySerialSymDenseMatrix.h"
%pythoncode
%{
class SerialSymDenseMatrix(UserArray,NumPySerialSymDenseMatrix):
    def __init__(self, *args):
      	"""
      	__init__(self) -> SerialSymDenseMatrix
      	__init__(self, PyObject array) -> SerialSymDenseMatrix
      	__init__(self, SerialSymDenseMatrix source) -> SerialSymDenseMatrix
      	"""
        NumPySerialSymDenseMatrix.__init__(self, *args)
        self.__initArray__()
    def __initArray__(self):
        self.array = self.A()
        self.__protected = True
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
        # This should get called when the SerialSymDenseMatrix is accessed after
        # not properly being initialized
        if not "array" in self.__dict__:
            self.__initArray__()
        try:
            return self.array.__getattribute__(key)
        except AttributeError:
            return SerialSymDenseMatrix.__getattribute__(self, key)
    def __setattr__(self, key, value):
        "Handle 'this' attribute properly and protect the 'array' and 'shape' attributes"
        if key == "this":
            NumPySerialSymDenseMatrix.__setattr__(self, key, value)
        else:
            if key in self.__dict__:
                if self.__protected:
                    if key == "array":
                        raise AttributeError, \
                              "Cannot change Epetra.SerialSymDenseMatrix array attribute"
                    if key == "shape":
                        raise AttributeError, \
                              "Cannot change Epetra.SerialSymDenseMatrix shape attribute"
            UserArray.__setattr__(self, key, value)
    def __getitem__(self, index):
        """
        __getitem__(self,int,int) -> int
        __getitem__(self,int,slice) -> array
        __getitem__(self,slice,int) -> array
        __getitem__(self,slice,slice) -> array
        """
        return self.array[index]
    def Shape(self,numRows,numCols):
        "Shape(self, int numRows, int numCols) -> int"
        result = NumPySerialSymDenseMatrix.Shape(self,numRows,numCols)
        self.__protected = False
        self.__initArray__()
        return result
    def Reshape(self,numRows,numCols):
        "Reshape(self, int numRows, int numCols) -> int"
        result = NumPySerialSymDenseMatrix.Reshape(self,numRows,numCols)
        self.__protected = False
        self.__initArray__()
        return result
_Epetra.NumPySerialSymDenseMatrix_swigregister(SerialSymDenseMatrix)
%}

//////////////////////////////////////
// Epetra_SerialDenseVector support //
//////////////////////////////////////
%ignore Epetra_SerialDenseVector::operator()(int);
%ignore Epetra_SerialDenseVector::operator()(int) const;
%inline
{
  struct SerialDenseVector{ };
}
%include "Epetra_SerialDenseVector.h"

///////////////////////////////////////////
// Epetra_NumPySerialDenseVector support //
///////////////////////////////////////////
%rename(NumPySerialDenseVector) Epetra_NumPySerialDenseVector;
%include "Epetra_NumPySerialDenseVector.h"
%pythoncode
%{
class SerialDenseVector(UserArray,NumPySerialDenseVector):
    def __init__(self, *args):
      	"""
      	__init__(self) -> SerialDenseVector
      	__init__(self, int length) -> SerialDenseVector
      	__init__(self, PyObject array) -> SerialDenseVector
      	__init__(self, SerialDenseVector source) -> SerialDenseVector
      	"""
        NumPySerialDenseVector.__init__(self, *args)
        self.__initArray__()
    def __initArray__(self):
        self.array = self.Values()
        self.__protected = True
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
        # This should get called when the SerialDenseVector is accessed after
        # not properly being initialized
        if not "array" in self.__dict__:
            self.__initArray__()
        try:
            return self.array.__getattribute__(key)
        except AttributeError:
            return SerialDenseVector.__getattribute__(self, key)
    def __setattr__(self, key, value):
        "Handle 'this' attribute properly and protect the 'array' attribute"
        if key == "this":
            NumPySerialDenseVector.__setattr__(self, key, value)
        else:
            if key in self.__dict__:
                if self.__protected:
                    if key == "array":
                        raise AttributeError, \
                              "Cannot change Epetra.SerialDenseVector array attribute"
            UserArray.__setattr__(self, key, value)
    def __call__(self,i):
        "__call__(self, int i) -> double"
        return self.__getitem__(i)
    def Size(self,length):
        "Size(self, int length) -> int"
        result = NumPySerialDenseVector.Size(self,length)
        self.__protected = False
        self.__initArray__()
        return result
    def Resize(self,length):
        "Resize(self, int length) -> int"
        result = NumPySerialDenseVector.Resize(self,length)
        self.__protected = False
        self.__initArray__()
        return result
_Epetra.NumPySerialDenseVector_swigregister(SerialDenseVector)
%}

//////////////////////////////////////
// Epetra_SerialDenseSolver support //
//////////////////////////////////////
%teuchos_rcp(Epetra_SerialDenseSolver)
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
%teuchos_rcp(Epetra_SerialDenseSVD)
%rename(SerialDenseSVD) Epetra_SerialDenseSVD;
%include "Epetra_SerialDenseSVD.h"

/////////////////////////////////////////
// Epetra_SerialSpdDenseSolver support //
/////////////////////////////////////////
// *** Epetra_SerialSpdDenseSolver is apparently not built ***
//#include "Epetra_SerialSpdDenseSolver.h"
//%teuchos_rcp(Epetra_SerialSpdDenseSolver)
//%rename(SerialSpdDenseSolver  ) Epetra_SerialSpdDenseSolver;
//%include "Epetra_SerialSpdDenseSolver.h"

%feature("compactdefaultargs");       // Turn the feature back on
