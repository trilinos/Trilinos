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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
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

      //! Create an empty matrix
      Matrix();

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

      //! Scale the matrix by a constant value.
      /*! This is needed to manipulate matrices in the LOCA library 
	routines. 
      */
      void scale(double value);

      //! Print this object to the given stream
      ostream& leftshift(ostream& stream) const;

      //! Prints out the matrix to the cout stream. 
      bool print() const;

      //! Returns the number of rows in the matrix
      int numRows() const;

      //! Returns the number of columns in the matrix
      int numCols() const;

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
