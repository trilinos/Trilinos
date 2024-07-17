// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
      Matrix(int m, int n) :
        // Compute size as size_t to avoid overflow for large matrices.
        p(m), q(n), entries(size_t(m)*size_t(n)) {}

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
        return ! stream.fail ();
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
