%{
// Epetra includes
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

// Local includes
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
%}

// Ignore directives
%ignore Epetra_MultiVector::operator=(const Epetra_MultiVector &);
%ignore Epetra_MultiVector::operator[](int);
%ignore Epetra_MultiVector::operator[](int) const;
%ignore Epetra_MultiVector::operator()(int) const;
%ignore Epetra_MultiVector::ExtractCopy(double *, int   ) const;  // These Extract methods
%ignore Epetra_MultiVector::ExtractCopy(double **       ) const;  // are given functionality
%ignore Epetra_MultiVector::ExtractView(double **, int *) const;  // in the derived class
%ignore Epetra_MultiVector::ExtractView(double ***      ) const;  // Epetra_NumPyMultiVector
%ignore Epetra_MultiVector::ResetView(double **);     // These are expert
%ignore Epetra_MultiVector::Pointers() const;         // method not supported in python
%ignore Epetra_Vector::ExtractCopy(double * ) const;  // These Extract methods are given functionality
%ignore Epetra_Vector::ExtractView(double **) const;  // in the derived class Epetra_NumPyMultiVector

// Rename directives
%rename(NumPyMultiVector) Epetra_NumPyMultiVector;
%rename(NumPyVector     ) Epetra_NumPyVector;

// Import directives for Epetra
%include "Epetra_MultiVector.h"
%include "Epetra_Vector.h"

// Local interface includes
%include "Epetra_NumPyMultiVector.h"
%include "Epetra_NumPyVector.h"

// Python code.  Here we define the Epetra.MultiVector and
// Epetra.Vector python classes, which multiply inherit from the
// Numeric UserArray class (making these classes Numeric arrays) and
// the Epetra_NumPyMultiVector or Epetra_NumPyVector class (making
// these classes also Epetra objects).
%pythoncode %{

from UserArray import *

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
        UserArray.__init__(self,self.ExtractView(),'d',copy=False,savespace=False)
    def __str__(self):
        return str(self.array)
    def __setattr__(self, key, value):
        "Protect the 'array' and 'shape' attributes"
        if key == "array":
            if key in self.__dict__:
                raise AttributeError, "Cannot change Epetra.MultiVector array attribute"
        elif key == "shape":
            value = tuple(value)
            if len(value) < 2:
                raise ValueError, "Epetra.MultiVector shape is " + str(value) + \
		  " but must have minimum of 2 elements"
        UserArray.__setattr__(self, key, value)

class Vector(UserArray,NumPyVector):
    def __init__(self, *args):
        NumPyVector.__init__(self, *args)
        UserArray.__init__(self,self.ExtractView(),'d',copy=False,savespace=False)
    def __str__(self):
        return str(self.array)
    def __setattr__(self, key, value):
        "Protect the 'array' attribute"
        if key == "array":
            if key in self.__dict__:
                raise AttributeError, "Cannot change Epetra.Vector array attribute"
        UserArray.__setattr__(self, key, value)

%}

%extend Epetra_MultiVector {
  bool isMultiVector() {
    return(true);
  }

  double __getitem__(int i, int j)
  {
    return((*self)[i][j]);
  }

  void __setitem__(int i, int j, double val)
  {
    (*self)[i][j] = val;
  }
}

%extend Epetra_Vector {
  bool isMultiVector() {
    return(false);
  }

  double __getitem__(int i)
  {
    return((*self)[i]);
  }

  void __setitem__(int i, double val)
  {
    (*self)[i] = val;
  }
}

