// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef NOX_EXAMPLE_MATRIX_H
#define NOX_EXAMPLE_MATRIX_H

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H" // for CopyType

namespace NOX {
namespace Example {

//! A simple square matrix class for use by NOX::Example::Group.
/*! The matrix is stored as a vector<double> array. */
class Matrix {

public:
  //! Create a square n x n matrix with all entries zero
  Matrix(int n);

  //! Copy constructor
  Matrix(const Matrix& a, CopyType type = NOX::DeepCopy);

  //! Destructor
  ~Matrix();

  //! Access the (i,j) entry of A
  double& operator()(int i, int j);

  //! Access the (i,j) entry of A
  const double& operator()(int i, int j) const;

private:

  //! This is a p x p matrix
  int p;

  //! Entries of the matrix
  vector<double> entries;
  
};

} // namespace Example
} // namespace NOX


#endif
