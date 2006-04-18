// @HEADER
// ***********************************************************************
//
//                 PyTrilinos: Rapid Prototyping Package
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

#ifndef NUMPYARRAY_H
#define NUMPYARRAY_H

#define NO_IMPORT_ARRAY
#include "numpy_include.h"

#include <iostream>

class NumPyArrayBase
{
public:
  NumPyArrayBase();
  virtual ~NumPyArrayBase () = 0;

  virtual std::ostream & print(std::ostream & stream) const = 0; 

  int                   getNumDims    () const;
  int                   getTotalLength() const;
  const int           * getDimLengths () const;
  double              * getDataArray  ()      ;
  const double        * getDataArray  () const;
  const int           * getStrides    () const;
  const PyArrayObject * getArrayObject() const;
  PyArrayObject       * getArrayObject()      ;
  
  // Reports if wrapped array is contiguous
  bool isContiguous() const;

protected:
  // Print base class info.  To be called by implementation
  // of print.
  std::ostream & printBase(std::ostream & stream) const; 

  void setData(PyArrayObject * p_pyArrayObject       ,
               int             numPyArrayNumDims     ,
               int             numPyArrayTotalLength ,
               int           * p_numPyArrayDimLengths,
               int           * p_numPyArrayStrides   ,
               double        * p_dataArray            );

private:
  // Private and not implemented
  NumPyArrayBase(const NumPyArrayBase & a_ref);
  const NumPyArrayBase & operator = (const NumPyArrayBase & a_rhs);

private:
  PyArrayObject * mp_pyArrayObject       ;
  int             m_numPyArrayNumDims    ;
  int             m_numPyArrayTotalLength;
  int           * mp_numPyArrayDimLengths;
  int           * mp_numPyArrayStrides   ;
  double        * mp_dataArray           ;
};

class NumPyArrayContiguous: public NumPyArrayBase
{
public:
  NumPyArrayContiguous(PyObject * p_pyObject);
  virtual ~NumPyArrayContiguous ();

  virtual std::ostream & print(std::ostream & stream) const; 
  
private:
  // Private and not implemented
  NumPyArrayContiguous();
  NumPyArrayContiguous(const NumPyArrayContiguous & a_ref);
  const NumPyArrayContiguous & operator = (const NumPyArrayContiguous & a_rhs);
};

class NumPyArray: public NumPyArrayBase
{
public:
  NumPyArray(PyObject * p_pyObject);
  virtual ~NumPyArray ();

  virtual std::ostream & print(std::ostream & stream) const; 

private:
  // Private and not implemented
  NumPyArray();
  NumPyArray(const NumPyArray & a_ref);
  const NumPyArray & operator = (const NumPyArray & a_rhs);

};

#endif // NumPyArray_h
