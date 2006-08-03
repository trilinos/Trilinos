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

#ifndef NOX_LAPACK_LINEARSOLVER_H
#define NOX_LAPACK_LINEARSOLVER_H

#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"

#include "NOX_LAPACK_Matrix.H"

namespace NOX {

  namespace LAPACK {

    //! A simple linear solver for use by NOX::LAPACK::Group
    /*!
     * This class provides a simple linear solver class that stores a 
     * NOX::LAPACK::Matrix and provides routines to apply the matrix
     * and solve it using BLAS and LAPACK routines.  It is templated
     * so that it can be used to solve both real and complex matrices.
     * It also stores an LU factorization of the matrix so repeated 
     * solves are more efficient.  The group should signal that the
     * matrix has changed by calling reset().
     */
    template <typename T>
    class LinearSolver {

    public:

      //! Constructor
      LinearSolver(int n);

      //! Copy constructor
      LinearSolver(const LinearSolver<T>& s);

      //! Destructor
      ~LinearSolver();

      //! Assignment
      LinearSolver& operator=(const LinearSolver<T>& s);

      //! Return matrix
      Matrix<T>& getMatrix();

      //! Return matrix
      const Matrix<T>& getMatrix() const;

      //! Reset
      /*!
       * Resets the LU factorization indicating the matrix is updated.
       */
      void reset();

      //! Apply matrix
      /*!
       * Set \c trans to \c true to apply the transpose.  \c ncols is the 
       * number of columns in \c input and \c output, which should be stored 
       * column-wise.
       */
      void apply(bool trans, int ncols, const T* input, T* output) const;

      //! Solve inverse
      /*!
       * Set \c trans to \c true to solve the transpose.  \c ncols is the 
       * number of columns in \c output, which should be stored column-
       * wise.  The right-hand-side is passed through \c output, which is
       * then overwritten with the result.
       */
      bool solve(bool trans, int ncols, T* output);

    protected:

      //! Matrix
      Matrix<T> mat;

      //! LU factorization of matrix
      Matrix<T> lu;

      //! Pivots
      std::vector<int> pivots;

      //! Is LU factorization valid
      bool isValidLU;

      //! BLAS wrappers
      Teuchos::BLAS<int,T> blas;

      //! LAPACK wrappers
      Teuchos::LAPACK<int,T> lapack;
      
    };

  } // namespace LAPACK

} // namespace NOX

template <typename T>
NOX::LAPACK::LinearSolver<T>::LinearSolver(int n) :
  mat(n,n),
  lu(n,n),
  pivots(n),
  isValidLU(false),
  blas(),
  lapack()
{
}

template <typename T>
NOX::LAPACK::LinearSolver<T>::LinearSolver(const NOX::LAPACK::LinearSolver<T>& s) :
  mat(s.mat),
  lu(s.lu),
  pivots(s.pivots),
  isValidLU(s.isValidLU),
  blas(),
  lapack()
{
}

template <typename T>
NOX::LAPACK::LinearSolver<T>::~LinearSolver()
{
}

template <typename T>
NOX::LAPACK::LinearSolver<T>&
NOX::LAPACK::LinearSolver<T>::operator=(const NOX::LAPACK::LinearSolver<T>& s)
{
  if (this != &s) {
    mat = s.mat;
    lu = s.lu;
    pivots = s.pivots;
    isValidLU = s.isValidLU;
  }

  return *this;
}

template <typename T>
NOX::LAPACK::Matrix<T>&
NOX::LAPACK::LinearSolver<T>::getMatrix()
{
  return mat;
}

template <typename T>
const NOX::LAPACK::Matrix<T>&
NOX::LAPACK::LinearSolver<T>::getMatrix() const
{
  return mat;
}

template <typename T>
void
NOX::LAPACK::LinearSolver<T>::reset()
{
  isValidLU = false;
}

template <typename T>
void
NOX::LAPACK::LinearSolver<T>::apply(bool trans, int ncols, const T* input, 
			      T* output) const
{
  Teuchos::ETransp tr = Teuchos::NO_TRANS;
  if (trans)
    tr = Teuchos::TRANS;

  int n = mat.numRows();
  blas.GEMM(tr, Teuchos::NO_TRANS, n, ncols, n, 1.0, &mat(0,0), n, 
	    input, n, 0.0, output, n);
}

template <typename T>
bool
NOX::LAPACK::LinearSolver<T>::solve(bool trans, int ncols, T* output)
{
  int info;
  int n = mat.numRows();

  // Compute LU factorization if necessary
  if (!isValidLU) {
    lu = mat;
    lapack.GETRF(n, n, &lu(0,0), n, &pivots[0], &info);
    if (info != 0)
      return false;
    isValidLU = true;
  }

  // Solve using LU factorization
  char tr = 'N';
  if (trans)
    tr = 'T';
  lapack.GETRS(tr, n, ncols, &lu(0,0), n, &pivots[0], output, n, &info);

  if (info != 0)
    return false;
  return true;
}

#endif
