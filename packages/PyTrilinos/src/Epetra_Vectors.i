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

// PyTrilinos configuration
%include "PyTrilinos_config.h"

///////////////////////////////////////////////////////////////////
// Handle the following Epetra_* and Epetra_NumPy* classes using //
// Teuchos::RCP<>                                                //
///////////////////////////////////////////////////////////////////
%teuchos_rcp_epetra_numpy(IntVector  )
%teuchos_rcp_epetra_numpy(MultiVector)
%teuchos_rcp_epetra_numpy(Vector     )
%teuchos_rcp_epetra_numpy(FEVector   )

//////////////////////////////
// Epetra_IntVector support //
//////////////////////////////
%inline
{
  struct IntVector{ };
}
%include "Epetra_IntVector.h"

////////////////////////////////
// Epetra_MultiVector support //
////////////////////////////////
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
%inline
{
  struct MultiVector{ };
}
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
%inline
{
  struct Vector{ };
}
%include "Epetra_Vector.h"

/////////////////////////////
// Epetra_FEVector support //
/////////////////////////////
%ignore Epetra_FEVector::ReplaceGlobalValues(int,int*,double*);
%ignore Epetra_FEVector::SumIntoGlobalValues(int,int*,double*);
%inline
{
  struct FEVector{ };
}
%include "Epetra_FEVector.h"

///////////////////////////////////
// Epetra_NumPyIntVector support //
///////////////////////////////////
namespace PyTrilinos
{
%feature("docstring")
Epetra_NumPyIntVector::ExtractCopy
"Return a numpy.ndarray that is a copy of the IntVector."
%feature("docstring")
Epetra_NumPyIntVector::ExtractView
"Return a numpy.ndarray that is a view of the IntVector."
%feature("docstring")
Epetra_NumPyIntVector::Values
"Return a numpy.ndarray that is a view of the IntVector."
%rename(NumPyIntVector) Epetra_NumPyIntVector;
// Wrappers for the Epetra_NumPyIntVector one- and two-argument
// constructors behave differently depending upon the version of swig
// (and possibly numpy) being used.  (Specifically, the generated
// dispatch functions assign priorities differently, and PyObject*
// arguments can have too high a priority.)  To avoid problems, I take
// control of these constructors here.
%extend Epetra_NumPyIntVector
{
  Epetra_NumPyIntVector(PyObject * arg1)
  {
    int                                 res   = 0;
    Epetra_IntVector                  * eiv   = NULL;
    PyTrilinos::Epetra_NumPyIntVector * enpiv = NULL;
    Epetra_BlockMap                   * bmap  = NULL;

#ifdef HAVE_TEUCHOS
    Teuchos::RCP< const Epetra_BlockMap > rcpbmap;
    void * vtemp  = NULL;
    int    newmem = 0;
    res = SWIG_ConvertPtrAndOwn(arg1, &vtemp, SWIGTYPE_p_Teuchos__RCPT_Epetra_BlockMap_t,
				0, &newmem);
    if (vtemp)
    {
      if (newmem & SWIG_CAST_NEW_MEMORY)
      {
	rcpbmap = *reinterpret_cast< Teuchos::RCP< const Epetra_BlockMap > * >(vtemp);
	delete reinterpret_cast< Teuchos::RCP< const Epetra_BlockMap > * >(vtemp);
	bmap = const_cast< Epetra_BlockMap * >(rcpbmap.get());
      }
      else
      {
	bmap = const_cast< Epetra_BlockMap * >(reinterpret_cast< Teuchos::RCP< const Epetra_BlockMap > * >(vtemp)->get());
      }
    }
#else
    res = SWIG_ConvertPtr(arg1, (void**)&bmap, SWIGTYPE_p_Epetra_BlockMap, 0);
#endif

    if (SWIG_CheckState(res))
    {
      enpiv = new PyTrilinos::Epetra_NumPyIntVector(*bmap);
    }
    else
    {
#ifdef HAVE_TEUCHOS
      Teuchos::RCP< const Epetra_IntVector > rcpeiv;
      vtemp = NULL;
      newmem = 0;
      res = SWIG_ConvertPtrAndOwn(arg1, &vtemp, SWIGTYPE_p_Teuchos__RCPT_Epetra_IntVector_t, 
				  0, &newmem);
      if (vtemp)
      {
	if (newmem & SWIG_CAST_NEW_MEMORY) 
	{
	  rcpeiv = *reinterpret_cast< Teuchos::RCP< const Epetra_IntVector > * >(vtemp);
	  delete reinterpret_cast< Teuchos::RCP< const Epetra_IntVector > * >(vtemp);
	  eiv = const_cast< Epetra_IntVector * >(rcpeiv.get());
	}
	else
	{
	  eiv = const_cast< Epetra_IntVector * >(reinterpret_cast< Teuchos::RCP< const Epetra_IntVector > * >(vtemp)->get());
        }
      }
#else
      res = SWIG_ConvertPtr(arg1, (void**)&eiv, SWIGTYPE_p_Epetra_IntVector, 0);
#endif

      if (SWIG_CheckState(res))
      {
	enpiv = new PyTrilinos::Epetra_NumPyIntVector(*eiv);
      }
      else
      {
	enpiv = new PyTrilinos::Epetra_NumPyIntVector(arg1);
      }
    }
    if (enpiv == NULL)
    {
      PyErr_SetString(PyExc_ValueError,
		      "Error constructing Epetra_NumPyIntVector\n"
		      "  Valid one-argument constructors:\n"
		      "    Epetra_NumPyIntVector(Epetra_BlockMap)\n"
		      "    Epetra_NumPyIntVector(Epetra_IntVector)\n"
		      "    Epetra_NumPyIntVector(array)");
    }
    return enpiv;
  }

  Epetra_NumPyIntVector(PyObject * arg1, PyObject * arg2)
  {
    int                                 res     = 0;
    bool                                zeroOut = true;
    PyTrilinos::Epetra_NumPyIntVector * enpiv   = NULL;
    Epetra_BlockMap                   * bmap    = NULL;

#ifdef HAVE_TEUCHOS
    Teuchos::RCP< const Epetra_BlockMap > rcpbmap;
    void * vtemp  = NULL;
    int    newmem = 0;
    res = SWIG_ConvertPtrAndOwn(arg1, &vtemp, SWIGTYPE_p_Teuchos__RCPT_Epetra_BlockMap_t,
				0, &newmem);
    if (vtemp)
    {
      if (newmem & SWIG_CAST_NEW_MEMORY)
      {
        rcpbmap = *reinterpret_cast< Teuchos::RCP< const Epetra_BlockMap > * >(vtemp);
        delete reinterpret_cast< Teuchos::RCP< const Epetra_BlockMap > * >(vtemp);
        bmap = const_cast< Epetra_BlockMap * >(rcpbmap.get());
      }
      else
      {
        bmap = const_cast< Epetra_BlockMap * >(reinterpret_cast< Teuchos::RCP< const Epetra_BlockMap > * >(vtemp)->get());
      }
    }
#else
    res = SWIG_ConvertPtr(arg1, (void**)&bmap, SWIGTYPE_p_Epetra_BlockMap, 0);
#endif

    if (SWIG_CheckState(res))
    {
      if (PyBool_Check(arg2))
      {
	zeroOut = (arg2 == Py_True) ? true : false;
	enpiv = new PyTrilinos::Epetra_NumPyIntVector(*bmap, zeroOut);
      }
      else
      {
	enpiv = new PyTrilinos::Epetra_NumPyIntVector(*bmap, arg2);
      }
    }
    if (enpiv == NULL)
    {
      PyErr_SetString(PyExc_ValueError,
		      "Error constructing Epetra_NumPyIntVector\n"
		      "  Valid two-argument constructors:\n"
		      "    Epetra_NumPyIntVector(Epetra_BlockMap, bool)\n"
		      "    Epetra_NumPyIntVector(Epetra_BlockMap, array)");
    }
    return enpiv;
  }
}
%ignore Epetra_NumPyIntVector::Epetra_NumPyIntVector(const Epetra_BlockMap&);
%ignore Epetra_NumPyIntVector::Epetra_NumPyIntVector(const Epetra_IntVector&);
%ignore Epetra_NumPyIntVector::Epetra_NumPyIntVector(PyObject*);
%ignore Epetra_NumPyIntVector::Epetra_NumPyIntVector(const Epetra_BlockMap&, bool);
%ignore Epetra_NumPyIntVector::Epetra_NumPyIntVector(const Epetra_BlockMap&, PyObject*);
}  // Namespace PyTrilinos
%include "Epetra_NumPyIntVector.h"
%pythoncode
%{
class IntVector(UserArray,NumPyIntVector):
    """
    Epetra.IntVector: A class for constructing and using dense integer vectors
    on a parallel computer.

    The Epetra.IntVector class enables the construction and use of integer dense
    vectors in a distributed memory environment. The distribution of the dense
    vector is determined in part by a Epetra.Comm object and a Epetra.Map (or
    Epetra.LocalMap or Epetra.BlockMap).

    Distributed Global vs. Replicated Local Distributed Global Vectors -
    In most instances, a multi-vector will be partitioned across multiple
    memory images associated with multiple processors. In this case, there
    is a unique copy of each element and elements are spread across all
    processors specified by the Epetra.Comm communicator.

    Replicated Local Vectors - Some algorithms use vectors that are too
    small to be distributed across all processors. Replicated local
    vectors handle these types of situation.

    In the python implementation, the IntVector stores an underlying numpy
    array, with which it shares the data buffer.  Also, almost all numpy array
    methods and operators are supported.
    """
    def __init__(self, *args):
        """
        __init__(self, BlockMap map, bool zeroOut=True) -> IntVector
        __init__(self, IntVector source) -> IntVector
        __init__(self, BlockMap map, PyObject array) -> IntVector
        __init__(self, PyObject array) -> IntVector

        Arguments:
            map      - BlockMap describing domain decomposition
            zeroOut  - Flag controlling whether to initialize IntVector to 0
            source   - Source IntVector for copy constructor
            array    - Python sequence that can be converted to a numpy array of
                       integers for initialization
        """
        NumPyIntVector.__init__(self, *args)
        self.__initArray__()
    def __initArray__(self):
        """
        __initArray__(self)
        
        Initialize the underlying numpy array.
        """
        UserArray.__init__(self, self.ExtractView(), dtype="i", copy=False)
    def __str__(self):
        """
        __str__(self)__ -> string
        
        Return a numpy-style string representation of the IntVector.
        """
        return str(self.array)
    def __lt__(self,other):
        """
        __lt__(self, other) -> bool

        Less-than operator (<).
        """
        return numpy.less(self.array,other)
    def __le__(self,other):
        """
        __le__(self, other) -> bool

        Less-than-or-equal operator (<=).
        """
        return numpy.less_equal(self.array,other)
    def __eq__(self,other):
        """
        __eq__(self, other) -> bool

        Equal operator (==).
        """
        return numpy.equal(self.array,other)
    def __ne__(self,other):
        """
        __ne__(self, other) -> bool

        Not-equal operator (!=).
        """
        return numpy.not_equal(self.array,other)
    def __gt__(self,other):
        """
        __gt__(self, other) -> bool

        Greater-than operator (>).
        """
        return numpy.greater(self.array,other)
    def __ge__(self,other):
        """
        __ge__(self, other) -> bool

        Greater-than-or-equal operator (>=).
        """
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
                    raise AttributeError, \
                          "Cannot change Epetra.IntVector array attribute"
            UserArray.__setattr__(self, key, value)
_Epetra.NumPyIntVector_swigregister(IntVector)
%}

/////////////////////////////////////
// Epetra_NumPyMultiVector support //
/////////////////////////////////////
namespace PyTrilinos
{
%feature("docstring")
Epetra_NumPyMultiVector::ExtractCopy
"Return a numpy.ndarray that is a copy of the MultiVector."
%feature("docstring")
Epetra_NumPyMultiVector::ExtractView
"Return a numpy.ndarray that is a view of the MultiVector."
%feature("docstring")
Epetra_NumPyMultiVector::Dot
"Return a numpy.ndarray of the dot products of the MultiVector and a."
%feature("docstring")
Epetra_NumPyMultiVector::Norm1
"Return a numpy.ndarray of the L-1 norms of MultiVector."
%feature("docstring")
Epetra_NumPyMultiVector::Norm2
"Return a numpy.ndarray of the the L-2 norms of MultiVector."
%feature("docstring")
Epetra_NumPyMultiVector::NormInf
"Return a numpy.ndarray of the L-infinity norms of MultiVector."
%feature("docstring")
Epetra_NumPyMultiVector::NormWeighted
"Return a numpy.ndarray of the weighted norms of MultiVector."
%feature("docstring")
Epetra_NumPyMultiVector::MinValue
"Return a numpy.ndarray of the minimum values in MultiVector."
%feature("docstring")
Epetra_NumPyMultiVector::MaxValue
"Return a numpy.ndarray of the maximum values in MultiVector."
%feature("docstring")
Epetra_NumPyMultiVector::MeanValue
"Return a numpy.ndarray of the mean values of the MultiVector."
%rename(NumPyMultiVector) Epetra_NumPyMultiVector;
}  // Namespace PyTrilinos
%include "Epetra_NumPyMultiVector.h"
%pythoncode
%{
class MultiVector(UserArray,NumPyMultiVector):
    """
    Epetra.MultiVector: A class for constructing and using dense multi- vectors,
    vectors and matrices in parallel.
    
    The Epetra.MultiVector class enables the construction and use of real-
    valued, double- precision dense vectors, multi-vectors, and matrices in a
    distributed memory environment. The dimensions and distribution of the dense
    multi-vectors is determined in part by a Epetra.Comm object, a Epetra.Map
    (or Epetra.LocalMap or Epetra.BlockMap) and the number of vectors passed to
    the constructors described below.
    
    There are several concepts that important for understanding the
    Epetra.MultiVector class:
    
    Multi-vectors, Vectors and Matrices. Vector - A list of real-valued,
    double-precision numbers. Also a multi-vector with one vector.
    
    Multi-Vector - A collection of one or more vectors, all having the same
    length and distribution.
    
    (Dense) Matrix - A special form of multi-vector such that stride in memory
    between any two consecutive vectors in the multi-vector is the same for all
    vectors. This is identical to a two-dimensional array in Fortran and plays
    an important part in high performance computations.
    
    Distributed Global vs. Replicated Local. Distributed Global Multi- vectors -
    In most instances, a multi-vector will be partitioned across multiple memory
    images associated with multiple processors. In this case, there is a unique
    copy of each element and elements are spread across all processors specified
    by the Epetra.Comm communicator.
    
    Replicated Local Multi-vectors - Some algorithms use multi-vectors that are
    too small to be distributed across all processors, the Hessenberg matrix in
    a GMRES computation. In other cases, such as with block iterative methods,
    block dot product functions produce small dense matrices that are required
    by all processors. Replicated local multi-vectors handle these types of
    situation.
    
    Multi-vector Functions vs. Dense Matrix Functions. Multi-vector functions -
    These functions operate simultaneously but independently on each vector in
    the multi-vector and produce individual results for each vector.
    
    Dense matrix functions - These functions operate on the multi-vector as a
    matrix, providing access to selected dense BLAS and LAPACK operations.

    In the python implementation, the MultiVector stores an underlying numpy
    array, with which it shares the data buffer.  This underlying numpy array
    has at least two dimensions, and the first dimension corresponds to the
    number of vectors.  Also, almost all numpy array methods and operators are
    supported.
    """
    def __init__(self, *args):
        """
        __init__(self, BlockMap map, int numVectors,
             bool zeroOut=True) -> MultiVector
        __init__(self, MultiVector source) -> MultiVector
        __init__(self, BlockMap map, PyObject array) -> MultiVector
        __init__(self, DataAccess CV, MultiVector source) -> MultiVector
        __init__(self, DataAccess CV, MultiVector source,
             PyObject range) -> MultiVector
        __init__(self, PyObject array) -> MultiVector

        Arguments:
            map         - BlockMap describing domain decomposition
            numVectors  - Number of vectors
            zeroOut     - Flag controlling whether to initialize MultiVector to
                          zero
            source      - Source MultiVector for copy constructors
            array       - Python sequence that can be converted to a numpy array
                          of doubles for initialization
            CV          - Epetra.Copy or Epetra.View
            range       - Python sequence specifying range of vector indexes
        """
        NumPyMultiVector.__init__(self, *args)
        self.__initArray__()
    def __initArray__(self):
        """
        __initArray__(self)

        Initialize the underlying numpy array.
        """
        UserArray.__init__(self, self.ExtractView(), dtype="d", copy=False)
    def __expand_index__(self, index):
        result = [slice(None, None, None)] * len(self.shape)
        if isinstance(index, tuple):
            for i in range(len(index)):
                result[i] = index[i]
        else:
            result[0] = index
        return tuple(result)
    def __getitem__(self,index):
        """
        x.__getitem__(y) <==> x[y]
        """
        result = UserArray.__getitem__(self,index)
        # If the result is an array (not a scalar), then we must take steps to
        # ensure that the resulting MultiVector has an accurate BlockMap
        if hasattr(result,"__len__"):
            # Obtain the new global IDs by getting a slice (based on index) from
            # an array of the old global IDs.  Use the new global IDs to build a
            # new BlockMap, upon which the new result will be based.
            index           = self.__expand_index__(index)
            newIndex        = index[1:]
            oldShape        = self.shape[1:]
            oldMap          = self.Map()
            gids            = oldMap.MyGlobalElements()
            gids.shape      = oldShape
            elemSizes       = oldMap.ElementSizeList()
            elemSizes.shape = oldShape
            newMap          = BlockMap(-1,
                                       gids[newIndex].ravel(),
                                       elemSizes[newIndex].ravel(),
                                       oldMap.IndexBase(),
                                       self.Comm())
            newShape        = result.shape
            if not (isinstance(index[0],slice) or hasattr(index[0],"__len__")):
                newShape = (1,) + newShape
            rarray       = result.array.ravel()
            rarray.shape = newShape
            result       = MultiVector(newMap, rarray)
        return result
    def __getslice__(self, i, j):
        """
        x.__getslice__(i,j) <==> x[i:j]
        """
        return self.__getitem__(slice(i,j))
    def __str__(self):
        """
        __str__(self) -> string

        Return the numpy-style string representation of the MultiVector.
        """
        return str(self.array)
    def __lt__(self,other):
        """
        __lt__(self, other) -> bool

        Less-than operator (<).
        """
        return numpy.less(self.array,other)
    def __le__(self,other):
        """
        __le__(self, other) -> bool

        Less-than-or-equal operator (<=).
        """
        return numpy.less_equal(self.array,other)
    def __eq__(self,other):
        """
        __eq__(self, other) -> bool

        Equal operator (==).
        """
        return numpy.equal(self.array,other)
    def __ne__(self,other):
        """
        __ne__(self, other) -> bool

        Not-equal operator (!=).
        """
        return numpy.not_equal(self.array,other)
    def __gt__(self,other):
        """
        __gt__(self, other) -> bool

        Greater-than operator (>).
        """
        return numpy.greater(self.array,other)
    def __ge__(self,other):
        """
        __ge__(self, other) -> bool

        Greater-than or equal operator (>=).
        """
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
                    raise AttributeError, \
                          "Cannot change Epetra.MultiVector array attribute"
            elif key == "shape":
                value = tuple(value)
                if len(value) < 2:
                    raise ValueError, "Epetra.MultiVector shape is " + \
                          str(value) + " but must have minimum of 2 elements"
            UserArray.__setattr__(self, key, value)
_Epetra.NumPyMultiVector_swigregister(MultiVector)
%}

////////////////////////////////
// Epetra_NumPyVector support //
////////////////////////////////
namespace PyTrilinos
{
%feature("docstring")
Epetra_NumPyVector::ExtractCopy
"Return a numpy.ndarray that is a copy of the Vector."
%feature("docstring")
Epetra_NumPyVector::ExtractView
"Return a numpy.ndarray that is a view of the Vector."
%feature("docstring")
Epetra_NumPyVector::Dot
"Return the dot product of the Vector and a."
%feature("docstring")
Epetra_NumPyVector::Norm1
"Return the L-1 norm of Vector."
%feature("docstring")
Epetra_NumPyVector::Norm2
"Return the the L-2 norm of Vector."
%feature("docstring")
Epetra_NumPyVector::NormInf
"Return the L-infinity norm of Vector."
%feature("docstring")
Epetra_NumPyVector::NormWeighted
"Return the weighted norm of Vector."
%feature("docstring")
Epetra_NumPyVector::MinValue
"Return the minimum values in Vector."
%feature("docstring")
Epetra_NumPyVector::MaxValue
"Return the maximum values in Vector."
%feature("docstring")
Epetra_NumPyVector::MeanValue
"Return the mean value of the Vector."
%feature("docstring")
Epetra_NumPyVector::ReplaceGlobalValues
"Replace global values at specified index (and offset)"
%feature("docstring")
Epetra_NumPyVector::ReplaceMyValues
"Replace local values at specified index (and offset)"
%feature("docstring")
Epetra_NumPyVector::SumIntoGlobalValues
"Sum into global values at specified indices (and offset)"
%feature("docstring")
Epetra_NumPyVector::SumIntoMyValues
"Sum into local values at specified indices (and offset)"
%rename(NumPyVector) Epetra_NumPyVector;
// Wrappers for the Epetra_NumPyVector one- and two-argument
// constructors behave differently depending upon the version of swig
// (and possibly numpy) being used.  (Specifically, the generated
// dispatch functions assign priorities differently, and PyObject*
// arguments can have too high a priority.)  To avoid problems, I take
// control of these constructors here.
%extend Epetra_NumPyVector
{
  Epetra_NumPyVector(PyObject * arg1)
  {
    int                              res  = 0;
    Epetra_Vector                  * ev   = NULL;
    PyTrilinos::Epetra_NumPyVector * enpv = NULL;
    Epetra_BlockMap                * bmap = NULL;

#ifdef HAVE_TEUCHOS
    Teuchos::RCP< const Epetra_BlockMap > rcpbmap;
    void * vtemp  = NULL;
    int    newmem = 0;
    res = SWIG_ConvertPtrAndOwn(arg1, &vtemp, SWIGTYPE_p_Teuchos__RCPT_Epetra_BlockMap_t,
				0, &newmem);
    if (vtemp)
    {
      if (newmem & SWIG_CAST_NEW_MEMORY)
      {
        rcpbmap = *reinterpret_cast< Teuchos::RCP< const Epetra_BlockMap > * >(vtemp);
        delete reinterpret_cast< Teuchos::RCP< const Epetra_BlockMap > * >(vtemp);
        bmap = const_cast< Epetra_BlockMap * >(rcpbmap.get());
      }
      else
      {
        bmap = const_cast< Epetra_BlockMap * >(reinterpret_cast< Teuchos::RCP< const Epetra_BlockMap > * >(vtemp)->get());
      }
    }
#else
    res = SWIG_ConvertPtr(arg1, (void**)&bmap, SWIGTYPE_p_Epetra_BlockMap, 0);
#endif

    if (SWIG_CheckState(res))
    {
      enpv = new PyTrilinos::Epetra_NumPyVector(*bmap);
    }
    else
    {
#ifdef HAVE_TEUCHOS
      Teuchos::RCP< const Epetra_Vector > rcpev;
      vtemp  = NULL;
      newmem = 0;
      res = SWIG_ConvertPtrAndOwn(arg1, &vtemp, SWIGTYPE_p_Teuchos__RCPT_Epetra_Vector_t,
				  0, &newmem);
      if (vtemp)
      {
	if (newmem & SWIG_CAST_NEW_MEMORY)
	{
	  rcpev = *reinterpret_cast< Teuchos::RCP< const Epetra_Vector > * >(vtemp);
	  delete reinterpret_cast< Teuchos::RCP< const Epetra_Vector > * >(vtemp);
	  ev = const_cast< Epetra_Vector * >(rcpev.get());
	}
	else
	{
	  ev = const_cast< Epetra_Vector * >(reinterpret_cast< Teuchos::RCP< const Epetra_Vector > * >(vtemp)->get());
	}
      }
#else
      res = SWIG_ConvertPtr(arg1, (void**)&ev, SWIGTYPE_p_Epetra_Vector, 0);
#endif
      if (SWIG_CheckState(res))
      {
	enpv = new PyTrilinos::Epetra_NumPyVector(*ev);
      }
      else
      {
	enpv = new PyTrilinos::Epetra_NumPyVector(arg1);
      }
    }
    if (enpv == NULL)
    {
      PyErr_SetString(PyExc_ValueError,
		      "Error constructing Epetra_NumPyVector\n"
		      "  Valid one-argument constructors:\n"
		      "    Epetra_NumPyVector(Epetra_BlockMap)\n"
		      "    Epetra_NumPyVector(Epetra_Vector)\n"
		      "    Epetra_NumPyVector(array)");
    }
    return enpv;
  }

  Epetra_NumPyVector(PyObject * arg1, PyObject * arg2)
  {
    int                              res     = 0;
    bool                             zeroOut = true;
    Epetra_BlockMap                * bmap    = NULL;
    Epetra_Vector                  * ev      = NULL;
    PyTrilinos::Epetra_NumPyVector * enpv    = NULL;
    Epetra_DataAccess                cv;

#ifdef HAVE_TEUCHOS
    Teuchos::RCP< const Epetra_BlockMap > rcpbmap;
    void * vtemp  = NULL;
    int    newmem = 0;
    res = SWIG_ConvertPtrAndOwn(arg1, &vtemp, SWIGTYPE_p_Teuchos__RCPT_Epetra_BlockMap_t,
				0, &newmem);
    if (vtemp)
    {
      if (newmem & SWIG_CAST_NEW_MEMORY)
      {
        rcpbmap = *reinterpret_cast< Teuchos::RCP< const Epetra_BlockMap > * >(vtemp);
        delete reinterpret_cast< Teuchos::RCP< const Epetra_BlockMap > * >(vtemp);
        bmap = const_cast< Epetra_BlockMap * >(rcpbmap.get());
      }
      else
      {
        bmap = const_cast< Epetra_BlockMap * >(reinterpret_cast< Teuchos::RCP< const Epetra_BlockMap > * >(vtemp)->get());
      }
    }
#else
    res = SWIG_ConvertPtr(arg1, (void**)&bmap, SWIGTYPE_p_Epetra_BlockMap, 0);
#endif

    if (SWIG_CheckState(res))
    {
      if (PyBool_Check(arg2))
      {
	zeroOut = (arg2 == Py_True) ? true : false;
	enpv = new PyTrilinos::Epetra_NumPyVector(*bmap, zeroOut);
      }
      else
      {
	enpv = new PyTrilinos::Epetra_NumPyVector(*bmap, arg2);
      }
    }
    else
    {
      if (PyInt_Check(arg1))
      {
	cv  = static_cast< Epetra_DataAccess >(PyInt_AsLong(arg1));
#ifdef HAVE_TEUCHOS
	Teuchos::RCP< const Epetra_Vector > rcpev;
	vtemp  = NULL;
	newmem = 0;
	res = SWIG_ConvertPtrAndOwn(arg2, &vtemp, SWIGTYPE_p_Teuchos__RCPT_Epetra_Vector_t,
				    0, &newmem);
	if (vtemp)
	{
	  if (newmem & SWIG_CAST_NEW_MEMORY)
	  {
	    rcpev = *reinterpret_cast< Teuchos::RCP< const Epetra_Vector > * >(vtemp);
	    delete reinterpret_cast< Teuchos::RCP< const Epetra_Vector > * >(vtemp);
	    ev = const_cast< Epetra_Vector * >(rcpev.get());
	  }
	  else
	  {
	    ev = const_cast< Epetra_Vector * >(reinterpret_cast< Teuchos::RCP< const Epetra_Vector > * >(vtemp)->get());
	  }
	}
#else
	res = SWIG_ConvertPtr(arg2, (void**)&ev, SWIGTYPE_p_Epetra_Vector, 0);
#endif
	if (SWIG_CheckState(res))
	{
	  enpv = new PyTrilinos::Epetra_NumPyVector(cv, *ev);
	}
      }
    }
    if (enpv == NULL)
    {
      PyErr_SetString(PyExc_ValueError,
		      "Error constructing Epetra_NumPyVector\n"
		      "  Valid two-argument constructors:\n"
		      "    Epetra_NumPyVector(Epetra_BlockMap, bool)\n"
		      "    Epetra_NumPyVector(Epetra_BlockMap, array)\n"
		      "    Epetra_NumPyVector(Epetra_DataAccess, Epetra_Vector)");
    }
    return enpv;
  }
}
%ignore Epetra_NumPyVector::Epetra_NumPyVector(const Epetra_BlockMap&);
%ignore Epetra_NumPyVector::Epetra_NumPyVector(const Epetra_Vector&);
%ignore Epetra_NumPyVector::Epetra_NumPyVector(PyObject*);
%ignore Epetra_NumPyVector::Epetra_NumPyVector(const Epetra_BlockMap&, bool);
%ignore Epetra_NumPyVector::Epetra_NumPyVector(const Epetra_BlockMap&, PyObject*);
%ignore Epetra_NumPyVector::Epetra_NumPyVector(Epetra_DataAccess, const Epetra_Vector&);
}
%include "Epetra_NumPyVector.h"
%pythoncode
%{
class Vector(UserArray,NumPyVector):
    """
    Epetra.Vector: A class for constructing and using dense vectors on a
    parallel computer.
    
    The Epetra.Vector class enables the construction and use of real- valued,
    double- precision dense vectors in a distributed memory environment. The
    distribution of the dense vector is determined in part by a Epetra.Comm
    object and a Epetra.Map (or Epetra.LocalMap or Epetra.BlockMap).
    
    This class is derived from the Epetra.MultiVector class. As such, it has
    full access to all of the functionality provided in the Epetra.MultiVector
    class.
    
    Distributed Global vs. Replicated Local Distributed Global Vectors - In most
    instances, a multi-vector will be partitioned across multiple memory images
    associated with multiple processors. In this case, there is a unique copy of
    each element and elements are spread across all processors specified by the
    Epetra.Comm communicator.
    
    Replicated Local Vectors - Some algorithms use vectors that are too small to
    be distributed across all processors. Replicated local vectors handle these
    types of situation.

    In the python implementation, the Vector stores an underlying numpy array,
    with which it shares the data buffer.  Also, almost all numpy array methods
    and operators are supported.
    """
    def __init__(self, *args):
        """
        __init__(self, BlockMap map, bool zeroOut=True) -> Vector
        __init__(self, Vector source) -> Vector
        __init__(self, BlockMap map, PyObject array) -> Vector
        __init__(self, DataAccess CV, Vector source) -> Vector
        __init__(self, DataAccess CV, MultiVector source, PyObject index) -> Vector
        __init__(self, PyObject array) -> Vector

        Arguments:
            map         - BlockMap describing domain decomposition
            zeroOut     - Flag controlling whether to initialize MultiVector to
                          zero
            source      - Source Vector or MultiVector for copy constructors
            array       - Python sequence that can be converted to a numpy array
                          of doubles for initialization
            CV          - Epetra.Copy or Epetra.View
            index       - MultiVector vector index for copy constructor
        """
        NumPyVector.__init__(self, *args)
        self.__initArray__()
    def __initArray__(self):
        """
        __initArray__(self)

        Initialize the underlying numpy array.
        """
        UserArray.__init__(self, self.ExtractView(), dtype="d", copy=False)
    def __getitem__(self,index):
        """
        x.__getitem__(y) <==> x[y]
        """
        result = UserArray.__getitem__(self,index)
        # If the result is an array (not a scalar) then we must take steps to
        # ensure that the resulting Vector has an accurate BlockMap
        if hasattr(result,"__len__"):
            # Obtain the new global IDs by getting a slice (based on index) from
            # an array of the old global IDs.  Use the new global IDs to build a
            # new BlockMap, upon which the new result will be based.
            oldMap          = self.Map()
            gids            = oldMap.MyGlobalElements()
            gids.shape      = self.shape
            elemSizes       = oldMap.ElementSizeList()
            elemSizes.shape = self.shape
            newMap          = BlockMap(-1,
                                       gids[index].ravel(),
                                       elemSizes[index].ravel(),
                                       oldMap.IndexBase(),
                                       self.Comm())
            newShape        = result.shape
            rarray          = result.array.ravel()
            rarray.shape    = newShape
            result          = Vector(newMap, rarray)
        return result
    def __getslice__(self, i, j):
        """
        x.__getslice__(i,j) <==> x[i:j]
        """
        return self.__getitem__(slice(i,j))
    def __str__(self):
        """
        __str__(self) -> string

        Return the numpy-style string representation of the MultiVector.
        """
        return str(self.array)
    def __lt__(self,other):
        """
        __lt__(self, other) -> bool

        Less-than operator (<).
        """
        return numpy.less(self.array,other)
    def __le__(self,other):
        """
        __le__(self, other) -> bool

        Less-than-or-equal operator (<=).
        """
        return numpy.less_equal(self.array,other)
    def __eq__(self,other):
        """
        __eq__(self, other) -> bool

        Equal operator (==).
        """
        return numpy.equal(self.array,other)
    def __ne__(self,other):
        """
        __ne__(self, other) -> bool

        Not-equal operator (!=).
        """
        return numpy.not_equal(self.array,other)
    def __gt__(self,other):
        """
        __gt__(self, other) -> bool

        Greater-than operator (>).
        """
        return numpy.greater(self.array,other)
    def __ge__(self,other):
        """
        __ge__(self, other) -> bool

        Greater-than or equal operator (>=).
        """
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
                    raise AttributeError, \
                          "Cannot change Epetra.Vector array attribute"
            UserArray.__setattr__(self, key, value)
_Epetra.NumPyVector_swigregister(Vector)
%}

////////////////////////////////
// Epetra_NumPyFEVector support //
////////////////////////////////
namespace PyTrilinos
{
%feature("docstring")
Epetra_NumPyFEVector::ExtractCopy
"Return a numpy.ndarray that is a copy of the FEVector."
%feature("docstring")
Epetra_NumPyFEVector::ExtractView
"Return a numpy.ndarray that is a view of the FEVector."
%feature("docstring")
Epetra_NumPyFEVector::Dot
"Return the dot product of the FEVector and a."
%feature("docstring")
Epetra_NumPyFEVector::Norm1
"Return the L-1 norm of FEVector."
%feature("docstring")
Epetra_NumPyFEVector::Norm2
"Return the the L-2 norm of FEVector."
%feature("docstring")
Epetra_NumPyFEVector::NormInf
"Return the L-infinity norm of FEVector."
%feature("docstring")
Epetra_NumPyFEVector::NormWeighted
"Return the weighted norm of FEVector."
%feature("docstring")
Epetra_NumPyFEVector::MinValue
"Return the minimum values in FEVector."
%feature("docstring")
Epetra_NumPyFEVector::MaxValue
"Return the maximum values in FEVector."
%feature("docstring")
Epetra_NumPyFEVector::MeanValue
"Return the mean value of the FEVector."
%feature("docstring")
Epetra_NumPyFEVector::ReplaceGlobalValues
"Replace global values at specified index (and offset)"
%feature("docstring")
Epetra_NumPyFEVector::SumIntoGlobalValues
"Sum into global values at specified indices (and offset)"
%rename(NumPyFEVector) Epetra_NumPyFEVector;
}  // Namespace PyTrilinos
%include "Epetra_NumPyFEVector.h"
%pythoncode
%{
class FEVector(UserArray,NumPyFEVector):
    """
    Epetra Finite-Element Vector. This class inherits Epetra.MultiVector and
    thus provides all Epetra.MultiVector functionality, with one restriction:
    currently an Epetra.FEVector only has 1 internal vector.
    
    The added functionality provided by Epetra.FEVector is the ability to
    perform finite-element style vector assembly. It accepts sub-vector
    contributions, such as those that would come from element-load vectors,
    etc., and these sub-vectors need not be wholly locally owned.  In other
    words, the user can assemble overlapping data (e.g., corresponding to shared
    finite-element nodes). When the user is finished assembling their vector
    data, they then call the method Epetra.FEVector::GlobalAssemble() which
    gathers the overlapping data (all non-local data that was input on each
    processor) into the data- distribution specified by the map that the
    Epetra.FEVector is constructed with.
    
    Note: At the current time (Sept 6, 2002) the methods in this implementation
    assume that there is only 1 point associated with each map element. This
    limitation will be removed in the near future.

    In the python implementation, the FEVector stores an underlying numpy array,
    with which it shares the data buffer.  Also, almost all numpy array methods
    and operators are supported.
    """
    def __init__(self, *args):
        """
        __init__(self, BlockMap map, bool zeroOut=True) -> FEVector
        __init__(self, FEVector source) -> FEVector
        __init__(self, BlockMap map, PyObject array) -> FEVector
        __init__(self, DataAccess CV, Vector source) -> FEVector
        __init__(self, DataAccess CV, MultiVector source, PyObject index) -> FEVector
        __init__(self, PyObject array) -> FEVector

        Arguments:
            map         - BlockMap describing domain decomposition
            zeroOut     - Flag controlling whether to initialize MultiVector to
                          zero
            source      - Source Vector or MultiVector for copy constructors
            array       - Python sequence that can be converted to a numpy array
                          of doubles for initialization
            CV          - Epetra.Copy or Epetra.View
            index       - MultiVector vector index for copy constructor
        """
        NumPyFEVector.__init__(self, *args)
        self.__initArray__()
    def __initArray__(self):
        """
        __initArray__(self)

        Initialize the underlying numpy array.
        """
        UserArray.__init__(self, self.ExtractView(), dtype="d", copy=False)
    def __str__(self):
        """
        __str__(self) -> string

        Return the numpy-style string representation of the MultiVector.
        """
        return str(self.array)
    def __lt__(self,other):
        """
        __lt__(self, other) -> bool

        Less-than operator (<).
        """
        return numpy.less(self.array,other)
    def __le__(self,other):
        """
        __le__(self, other) -> bool

        Less-than-or-equal operator (<=).
        """
        return numpy.less_equal(self.array,other)
    def __eq__(self,other):
        """
        __eq__(self, other) -> bool

        Equal operator (==).
        """
        return numpy.equal(self.array,other)
    def __ne__(self,other):
        """
        __ne__(self, other) -> bool

        Not-equal operator (!=).
        """
        return numpy.not_equal(self.array,other)
    def __gt__(self,other):
        """
        __gt__(self, other) -> bool

        Greater-than operator (>).
        """
        return numpy.greater(self.array,other)
    def __ge__(self,other):
        """
        __ge__(self, other) -> bool

        Greater-than or equal operator (>=).
        """
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
