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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_LAPACK_MATRIX_H
#define NOX_LAPACK_MATRIX_H

#include "NOX_Common.H"

namespace NOX {

  namespace LAPACK {

    //! A simple square matrix class for use by NOX::LAPACK::Group.
    /*! 
     * The matrix is stored as a std::vector<T> array. It is templated so
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
std::ostream& operator<<(std::ostream& stream, const NOX::LAPACK::Matrix<T>& v) {
  v.print(stream);
  return stream;
}


#endif
