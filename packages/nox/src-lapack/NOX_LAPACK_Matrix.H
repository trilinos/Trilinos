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

#ifndef NOX_LAPACK_MATRIX_H
#define NOX_LAPACK_MATRIX_H

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H" // for CopyType
#include "NOX_LAPACK_Vector.H"

namespace NOX {
  namespace LAPACK {

    //! A simple square matrix class for use by NOX::LAPACK::Group.
    /*! The matrix is stored as a vector<double> array. */
    class Matrix {

    public:
      //! Create a m x n matrix with all entries zero
      Matrix(int m, int n);

      //! Copy constructor
      Matrix(const Matrix& a, CopyType type = NOX::DeepCopy);

      //! Destructor
      ~Matrix();

      //! Access the (i,j) entry of A
      double& operator()(int i, int j);

      //! Access the (i,j) entry of A
      const double& operator()(int i, int j) const;

      ostream& leftshift(ostream& stream) const;

      //! Prints out the matrix to the cout stream. 
      bool print() const;

    private:

      //! This is a p x q matrix
      int p, q;

      //! Entries of the matrix
      Vector entries;
  
    };

  } // namespace LAPACK
} // namespace NOX

ostream& operator<<(ostream& stream, const NOX::LAPACK::Matrix& v);


#endif
