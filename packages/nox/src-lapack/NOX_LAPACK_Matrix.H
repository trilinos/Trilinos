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

namespace NOX {

  namespace LAPACK {

    //! A simple square matrix class for use by NOX::LAPACK::Group.
    /*! 
     * The matrix is stored as a vector<T> array. It is templated so
     * it can store entries of different types.  For example, the LOCA
     * LAPACK group stores a complex matrix for Hopf tracking.
     */
    template <typename T>
    class Matrix {

    public:

      //! Create an empty matrix
      Matrix() : p(0), q(0), entries() {}

      //! Create a m x n matrix with all entries zero
      Matrix(int m, int n) : p(m), q(n), entries(m*n) {}

      //! Copy constructor
      Matrix(const Matrix& a) :
	p(a.p), q(a.q), entries(a.entries) {}

      //! Destructor
      ~Matrix() {}

      //! Access the (i,j) entry of A
      T& operator()(int i, int j) { return entries[i + (p*j)]; }

      //! Access the (i,j) entry of A
      const T& operator()(int i, int j) const { return entries[i + (p*j)]; }

      //! Scale the matrix by a constant value.
      /*! This is needed to manipulate matrices in the LOCA library 
	routines. 
      */
      void scale(T v ) {
	for (int i=0; i<p*q; i++)
	  entries[i] *= v;
      }

      //! Prints out the matrix 
      bool print(std::ostream& stream) const {
	for (int i=0; i<p; i++) {
	  stream << "[ ";
	  for (int j=0; j<q; j++)
	    stream << operator()(i,j) << " ";
	  stream << "]" << std::endl;
	}
	return stream;
      }

      //! Returns the number of rows in the matrix
      int numRows() const { return p; }

      //! Returns the number of columns in the matrix
      int numCols() const { return q; }

    private:

      //! This is a p x q matrix
      int p, q;

      //! Entries of the matrix
      std::vector<T> entries;
  
    };

  } // namespace LAPACK

} // namespace NOX

template <typename T>
ostream& operator<<(ostream& stream, const NOX::LAPACK::Matrix<T>& v) {
  v.print(stream);
  return stream;
}


#endif
