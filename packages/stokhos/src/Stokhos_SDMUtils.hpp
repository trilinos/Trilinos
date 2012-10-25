// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_SDM_UTILS_HPP
#define STOKHOS_SDM_UTILS_HPP

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_LAPACK.hpp"
#include <ostream>

#define DGEQPF_F77  F77_BLAS_MANGLE(dgeqpf,DGEQPF)
#define DGEQP3_F77  F77_BLAS_MANGLE(dgeqp3,DGEQP3)
extern "C" {
void DGEQPF_F77(int*, int*, double*, int*, int*, double*, double*, int*);
void DGEQP3_F77(int*, int*, double*, int*, int*, double*, double*, int*, int*);
}

#include "Stokhos_ConfigDefs.h"

#ifdef HAVE_STOKHOS_MATLABLIB
extern "C" {
#include "mat.h"
#include "matrix.h"
}
#endif

namespace Stokhos {

  namespace detail {

    //! Compute weighted inner product between vectors x and y
    template <typename ordinal_type, typename scalar_type>
    scalar_type weighted_inner_product(
      const Teuchos::SerialDenseVector<ordinal_type,scalar_type>& x,
      const Teuchos::SerialDenseVector<ordinal_type,scalar_type>& y,
      const Teuchos::Array<scalar_type>& w)
    {
      ordinal_type n = x.length();
      scalar_type t = 0;
      for (ordinal_type i=0; i<n; i++)
	t += x[i]*w[i]*y[i];
      return t;
    }

    //! Overwrite x with alpha*x + beta*y
    template <typename ordinal_type, typename scalar_type>
    void saxpy(
      const scalar_type& alpha,
      Teuchos::SerialDenseVector<ordinal_type, scalar_type>& x,
      const scalar_type& beta,
      const Teuchos::SerialDenseVector<ordinal_type, scalar_type>& y)
    {
      ordinal_type n = x.length();
      for (ordinal_type i=0; i<n; i++)
	x[i] = alpha*x[i] + beta*y[i];
    }

  }

  // Print a matrix in matlab-esque format
  template <typename ordinal_type, typename scalar_type>
  void
  print_matlab(std::ostream& os, 
	       const Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>& A)
  {
    os << "[ ";
    for (ordinal_type i=0; i<A.numRows(); i++) {
      for (ordinal_type j=0; j<A.numCols(); j++)
	os << A(i,j) << " ";
      if (i < A.numRows()-1)
	os << ";" << std::endl << "  ";
    }
    os << "];" << std::endl;  
  }

#ifdef HAVE_STOKHOS_MATLABLIB
  // Save a matrix to matlab MAT file
  template <typename ordinal_type>
  void
  save_matlab(const std::string& file_name, const std::string& matrix_name, 
	      const Teuchos::SerialDenseMatrix<ordinal_type,double>& A,
	      bool overwrite = true)
  {
    // Open matfile
    const char *mode;
    if (overwrite)
      mode = "w";
    else
      mode = "u";
    MATFile *file = matOpen(file_name.c_str(), mode);
    TEUCHOS_ASSERT(file != NULL);

    // Convert SDM to Matlab matrix
    mxArray *mat = mxCreateDoubleMatrix(0, 0, mxREAL);
    TEUCHOS_ASSERT(mat != NULL);
    mxSetPr(mat, A.values());
    mxSetM(mat, A.numRows());
    mxSetN(mat, A.numCols());

    // Save matrix to file
    ordinal_type ret = matPutVariable(file, matrix_name.c_str(), mat);
    TEUCHOS_ASSERT(ret == 0);

    // Close file
    ret = matClose(file);
    TEUCHOS_ASSERT(ret == 0);

    // Free matlab matrix
    mxSetPr(mat, NULL);
    mxSetM(mat, 0);
    mxSetN(mat, 0);
    mxDestroyArray(mat);
  }
#endif

  //! Vector-infinity norm of a matrix
  template <typename ordinal_type, typename scalar_type>
  scalar_type vec_norm_inf(
    const Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& A) {
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type> vec_A(
      Teuchos::View, A.values(), 1, A.numRows()*A.numCols(), 1);
    return vec_A.normInf();
  }

  //! Compute thin QR using classical Gram-Schmidt
  /*!
   * For A an m-by-n matrix computes A = Q*R with R k-by-k upper triangular,
   * Q m-by-k with orthogonal columns, k <= min(m,n).
   */
  template <typename ordinal_type, typename scalar_type>
  void QR_CGS(
    ordinal_type k,
    const Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>& A,
    const Teuchos::Array<scalar_type>& w,
    Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>& Q,
    Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>& R)
  {
    using Teuchos::getCol;
    typedef Teuchos::SerialDenseVector<ordinal_type,scalar_type> SDV;
    typedef Teuchos::SerialDenseMatrix<ordinal_type,scalar_type> SDM;

    // getCol() requires non-const SDM
    SDM& Anc = const_cast<SDM&>(A);

    // Make sure Q is the right size
    ordinal_type m = A.numRows();
    if (Q.numRows() != m || Q.numCols() != k)
      Q.shape(m,k);
    if (R.numRows() != k || R.numCols() != k)
      R.shape(k,k);
  
    for (ordinal_type j=0; j<k; j++) {
      SDV Aj = getCol(Teuchos::View, Anc, j);
      SDV Qj = getCol(Teuchos::View, Q, j);
      Qj.assign(Aj);
      for (ordinal_type i=0; i<j; i++) {
	SDV Qi = getCol(Teuchos::View, Q, i);
	R(i,j) = detail::weighted_inner_product(Qi, Aj, w);
	detail::saxpy(1.0, Qj, -R(i,j), Qi);  // Q(:,j) = 1.0*Q(:,j) - R(i,j)*Q(:,i)
      }
      R(j,j) = std::sqrt(detail::weighted_inner_product(Qj, Qj, w));
      Qj.scale(1.0/R(j,j));
    }

  }

  //! Compute thin QR using modified Gram-Schmidt
  /*!
   * For A an m-by-n matrix computes A = Q*R with R k-by-k upper triangular,
   * Q m-by-k with orthogonal columns, k <= min(m,n).
   */
  template <typename ordinal_type, typename scalar_type>
  void QR_MGS(
    ordinal_type k,
    const Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>& A,
    const Teuchos::Array<scalar_type>& w,
    Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>& Q,
    Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>& R)
  {
    using Teuchos::getCol;
    typedef Teuchos::SerialDenseVector<ordinal_type,scalar_type> SDV;
    typedef Teuchos::SerialDenseMatrix<ordinal_type,scalar_type> SDM;

    // getCol() requires non-const SDM
    SDM& Anc = const_cast<SDM&>(A);

    // Make sure Q is the right size
    ordinal_type m = A.numRows();
    if (Q.numRows() != m || Q.numCols() != k)
      Q.shape(m,k);
    if (R.numRows() != k || R.numCols() != k)
      R.shape(k,k);
    
    for (ordinal_type i=0; i<k; i++) {
      SDV Ai = getCol(Teuchos::View, Anc, i);
      SDV Qi = getCol(Teuchos::View, Q, i);
      Qi.assign(Ai);
    }
    for (ordinal_type i=0; i<k; i++) {
      SDV Qi = getCol(Teuchos::View, Q, i);
      for (ordinal_type j=0; j<i; j++) {
	SDV Qj = getCol(Teuchos::View, Q, j);
	scalar_type v = detail::weighted_inner_product(Qi, Qj, w);
	detail::saxpy(1.0, Qi, -v, Qj);  // Q(:,i) = 1.0*Q(:,i) - v*Q(:,j)
	R(j,i) += v;
      }
      R(i,i) = std::sqrt(detail::weighted_inner_product(Qi, Qi, w));
      Qi.scale(1.0/R(i,i));
    }
  }

  //! Compute thin QR using modified Gram-Schmidt with reorthogonalization
  /*!
   * For A an m-by-n matrix computes A = Q*R with R k-by-k upper triangular,
   * Q m-by-k with orthogonal columns, k <= min(m,n).
   */
  template <typename ordinal_type, typename scalar_type>
  void QR_MGS2(
    ordinal_type k,
    const Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>& A,
    const Teuchos::Array<scalar_type>& w,
    Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>& Q,
    Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>& R)
  {
    using Teuchos::getCol;
    typedef Teuchos::SerialDenseVector<ordinal_type,scalar_type> SDV;
    typedef Teuchos::SerialDenseMatrix<ordinal_type,scalar_type> SDM;

    // getCol() requires non-const SDM
    SDM& Anc = const_cast<SDM&>(A);

    // Make sure Q is the right size
    ordinal_type m = A.numRows();
    if (Q.numRows() != m || Q.numCols() != k)
      Q.shape(m,k);
    if (R.numRows() != k || R.numCols() != k)
      R.shape(k,k);
    
    for (ordinal_type i=0; i<k; i++) {
      SDV Ai = getCol(Teuchos::View, Anc, i);
      SDV Qi = getCol(Teuchos::View, Q, i);
      Qi.assign(Ai);
    }
    for (ordinal_type i=0; i<k; i++) {
      SDV Qi = getCol(Teuchos::View, Q, i);
      for (ordinal_type j=0; j<i; j++) {
	SDV Qj = getCol(Teuchos::View, Q, j);
	scalar_type v = detail::weighted_inner_product(Qi, Qj, w);
	detail::saxpy(1.0, Qi, -v, Qj);  // Q(:,i) = 1.0*Q(:,i) - v*Q(:,j)
	R(j,i) += v;
      }
      for (ordinal_type j=i-1; j>=0; j--) {
	SDV Qj = getCol(Teuchos::View, Q, j);
	scalar_type v = detail::weighted_inner_product(Qi, Qj, w);
	detail::saxpy(1.0, Qi, -v, Qj);  // Q(:,i) = 1.0*Q(:,i) - v*Q(:,j)
	R(j,i) += v;
      }
      R(i,i) = std::sqrt(detail::weighted_inner_product(Qi, Qi, w));
      Qi.scale(1.0/R(i,i));
    }
  }

  //! Compute thin QR using Householder reflections
  /*!
   * For A an m-by-n matrix computes A = Q*R with R k-by-k upper triangular,
   * Q m-by-k with orthogonal columns, k <= min(m,n).
   *
   * The QR factorization is computed by the corresponding LAPACK function.
   */
  template <typename ordinal_type, typename scalar_type>
  void
  QR_Householder(
    ordinal_type k,
    const Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& A,
    const Teuchos::Array<scalar_type>& w,
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& Q,
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& R) 
  {
    Teuchos::LAPACK<ordinal_type,scalar_type> lapack;
    ordinal_type m = A.numRows();
    ordinal_type n = A.numCols();
    ordinal_type kk = std::min(m,n);
    if (k > kk)
      k = kk;

    // Check that each component of w is 1
    for (ordinal_type i=0; i<w.size(); i++)
      TEUCHOS_TEST_FOR_EXCEPTION(
	w[i] != 1.0, std::logic_error, 
	"CPQR_Householder_threshold() requires unit weight vector!");

    // Lapack routine overwrites A, so copy into temporary matrix
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type> AA(
      Teuchos::Copy, A, m, n);

    // QR
    ordinal_type lda = AA.stride();
    Teuchos::Array<scalar_type> tau(k);
    Teuchos::Array<scalar_type> work(1);
    ordinal_type info;
    ordinal_type lwork = -1;
    lapack.GEQRF(m, n, AA.values(), lda, &tau[0], &work[0], lwork, &info);
    TEUCHOS_TEST_FOR_EXCEPTION(
      info < 0, std::logic_error, "geqrf returned info = " << info);
    lwork = work[0];
    work.resize(lwork);
    lapack.GEQRF(m, n, AA.values(), lda, &tau[0], &work[0], lwork, &info);
    TEUCHOS_TEST_FOR_EXCEPTION(
      info < 0, std::logic_error, "geqrf returned info = " << info);

    // Extract R
    if (R.numRows() != k || R.numCols() != n)
      R.shape(k,n);
    R.putScalar(0.0);
    for (ordinal_type i=0; i<k; i++)
      for (ordinal_type j=i; j<n; j++)
	R(i,j) = AA(i,j);

    // Extract Q
    if (Q.numRows() != m || Q.numCols() != k)
      Q.shape(m,k);
    lwork = -1;
    lapack.ORGQR(m, k, k, AA.values(), lda, &tau[0], &work[0], lwork, &info);
    TEUCHOS_TEST_FOR_EXCEPTION(
      info < 0, std::logic_error, "orgqr returned info = " << info);
    lwork = work[0];
    work.resize(lwork);
    lapack.ORGQR(m, k, k, AA.values(), lda, &tau[0], &work[0], lwork, &info);
    TEUCHOS_TEST_FOR_EXCEPTION(
      info < 0, std::logic_error, "orgqr returned info = " << info);
    if (Q.numRows() != m || Q.numCols() != k)
      Q.shape(m, k);
    for (ordinal_type i=0; i<m; i++)
      for (ordinal_type j=0; j<k; j++)
	Q(i,j) = AA(i,j);
  }

  //! Compute column-pivoted QR using Householder reflections
  /*!
   * For A an m-by-n matrix with m >= n, computes A*P = Q*R with R
   * n-by-n upper triangular and Q m-by-n with orthogonal columns (often
   * called the economy size QR) and P an m-by-n permutation matrix.  For
   * n >= m, computes A*P = Q*R with R m-by-n upper trapezoidal and Q
   * m-by-m upper trapezoidal (R = [R_1 R_2] with R_1 upper triangular and
   * R_2 rectangular).  For k = min(m,n), both cases are handled with 
   * Q m-by-k and R k-by-n.
   *
   * The QR factorization is computed by the corresponding LAPACK function.
   */
  template <typename ordinal_type, typename scalar_type>
  void
  CPQR_Householder(
    const Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& A,
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& Q,
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& R,
    Teuchos::Array<ordinal_type>& piv) 
  {
    Teuchos::LAPACK<ordinal_type,scalar_type> lapack;
    ordinal_type m = A.numRows();
    ordinal_type n = A.numCols();
    ordinal_type k = std::min(m,n);

    // Lapack routine overwrites A, so copy into temporary matrix
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type> AA(
      Teuchos::Copy, A, m, n);
    if (Q.numRows() != m || Q.numCols() != k)
      Q.shape(m,k);

    // Teuchos LAPACK interface doesn't have dgeqpf, so call it directly

    // Column pivoted QR
    ordinal_type lda = AA.stride();
    piv.resize(n);
    Teuchos::Array<scalar_type> tau(k);
    Teuchos::Array<scalar_type> work(3*n);
    ordinal_type info;
    DGEQPF_F77(&m, &n, AA.values(), &lda, &piv[0], &tau[0], &work[0], &info);
    TEUCHOS_TEST_FOR_EXCEPTION(
      info < 0, std::logic_error, "dgeqp3 returned info = " << info);

    // Extract R
    if (R.numRows() != k || R.numCols() != n)
      R.shape(k,n);
    R.putScalar(0.0);
    for (ordinal_type i=0; i<k; i++)
      for (ordinal_type j=i; j<n; j++)
	R(i,j) = AA(i,j);

    // Extract Q
    ordinal_type lwork = -1;
    lapack.ORGQR(m, k, k, AA.values(), lda, &tau[0], &work[0], lwork, &info);
    TEUCHOS_TEST_FOR_EXCEPTION(
      info < 0, std::logic_error, "orgqr returned info = " << info);
    lwork = work[0];
    work.resize(lwork);
    lapack.ORGQR(m, k, k, AA.values(), lda, &tau[0], &work[0], lwork, &info);
    TEUCHOS_TEST_FOR_EXCEPTION(
      info < 0, std::logic_error, "orgqr returned info = " << info);
    if (Q.numRows() != m || Q.numCols() != k)
      Q.shape(m, k);
    for (ordinal_type i=0; i<m; i++)
      for (ordinal_type j=0; j<k; j++)
	Q(i,j) = AA(i,j);

    // Transform piv to zero-based indexing
    for (ordinal_type i=0; i<n; i++)
      piv[i] -= 1;
  }

  //! Compute column-pivoted QR using Householder reflections
  /*!
   * For A an m-by-n matrix with m >= n, computes A*P = Q*R with R
   * n-by-n upper triangular and Q m-by-n with orthogonal columns (often
   * called the economy size QR) and P an m-by-n permutation matrix.  For
   * n >= m, computes A*P = Q*R with R m-by-n upper trapezoidal and Q
   * m-by-m upper trapezoidal (R = [R_1 R_2] with R_1 upper triangular and
   * R_2 rectangular).  For k = min(m,n), both cases are handled with 
   * Q m-by-k and R k-by-n.
   *
   * The QR factorization is computed by the corresponding LAPACK function.
   * This version uses the BLAS3-rich xGEQP3.
   */
  template <typename ordinal_type, typename scalar_type>
  void
  CPQR_Householder3(
    const Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& A,
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& Q,
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& R,
    Teuchos::Array<ordinal_type>& piv) 
  {
    Teuchos::LAPACK<ordinal_type,scalar_type> lapack;
    ordinal_type m = A.numRows();
    ordinal_type n = A.numCols();
    ordinal_type k = std::min(m,n);

    // Lapack routine overwrites A, so copy into temporary matrix
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type> AA(
      Teuchos::Copy, A, m, n);
    if (Q.numRows() != m || Q.numCols() != k)
      Q.shape(m,k);

    // Teuchos LAPACK interface doesn't have dgeqp3, so call it directly

    // Workspace query
    ordinal_type lda = AA.stride();
    piv.resize(n);
    Teuchos::Array<scalar_type> tau(k);
    Teuchos::Array<scalar_type> work(1);
    ordinal_type lwork = -1;
    ordinal_type info;
    DGEQP3_F77(&m, &n, AA.values(), &lda, &piv[0], &tau[0], &work[0], &lwork, 
	     &info);
    TEUCHOS_TEST_FOR_EXCEPTION(
      info < 0, std::logic_error, "dgeqp3 returned info = " << info);

    // Column pivoted QR
    lwork = work[0];
    work.resize(lwork);
    DGEQP3_F77(&m, &n, AA.values(), &lda, &piv[0], &tau[0], &work[0], &lwork, 
	       &info);
    TEUCHOS_TEST_FOR_EXCEPTION(
      info < 0, std::logic_error, "dgeqp3 returned info = " << info);

    // Extract R
    if (R.numRows() != k || R.numCols() != n)
      R.shape(k,n);
    R.putScalar(0.0);
    for (ordinal_type i=0; i<k; i++)
      for (ordinal_type j=i; j<n; j++)
	R(i,j) = AA(i,j);

    // Extract Q
    lwork = -1;
    lapack.ORGQR(m, k, k, AA.values(), lda, &tau[0], &work[0], lwork, &info);
    TEUCHOS_TEST_FOR_EXCEPTION(
      info < 0, std::logic_error, "orgqr returned info = " << info);
    lwork = work[0];
    work.resize(lwork);
    lapack.ORGQR(m, k, k, AA.values(), lda, &tau[0], &work[0], lwork, &info);
    TEUCHOS_TEST_FOR_EXCEPTION(
      info < 0, std::logic_error, "orgqr returned info = " << info);
    if (Q.numRows() != m || Q.numCols() != k)
      Q.shape(m, k);
    for (ordinal_type i=0; i<m; i++)
      for (ordinal_type j=0; j<k; j++)
	Q(i,j) = AA(i,j);

    // Transform piv to zero-based indexing
    for (ordinal_type i=0; i<n; i++)
      piv[i] -= 1;
  }

  //! Compute column-pivoted QR using Householder reflections
  /*!
   * For A an m-by-n matrix, computes A*P = Q*R with R k-by-k upper triangular,
   * Q m-by-k with orthonormal columns, and P an n-by-k permutation matrix.
   * Here k <= min(m,n) is determined by a rank threshold tau provided by the
   * user.  The resulting R will have cond(R) <= 1/tau.  P is returned in the
   * pivot array \c piv and the rank k returned by the function.  Only the first
   * k entries of \c piv will be set.  As with LAPACK, the user can require
   * columns of A to be included in P by setting the corresponding entries
   * of piv to be nonzero on input.
   *
   * If \c make_R_square is \c false then R is k-by-n.
   *
   * This ultimately uses the LAPACK column-pivoted QR function which
   * does a full QR factorization.  This then extracts the parts of Q, R, and P
   * determined by the threshold as described above.  As such, this function
   * requires the weight vector to be 1 (Note the weight vector will be ignored
   * if it is size 0).
   */
  template <typename ordinal_type, typename scalar_type>
  ordinal_type
  CPQR_Householder_threshold(
    const scalar_type& rank_threshold,
    const Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& A,
    const Teuchos::Array<scalar_type>& w,
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& Q,
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& R,
    Teuchos::Array<ordinal_type>& piv) 
  {
    // Check that each component of w is 1
    for (ordinal_type i=0; i<w.size(); i++)
      TEUCHOS_TEST_FOR_EXCEPTION(
	w[i] != 1.0, std::logic_error, 
	"CPQR_Householder_threshold() requires unit weight vector!");

    // Compute full QR
    CPQR_Householder(A, Q, R, piv);

    // Find leading block of R such that cond(R) <= 1/tau
    ordinal_type rank = 0;
    ordinal_type m = R.numRows();
    scalar_type r_max = std::abs(R(rank,rank));
    scalar_type r_min = std::abs(R(rank,rank));
    for (rank=1; rank<m; rank++) {
      if (std::abs(R(rank,rank)) > r_max)
	r_max = std::abs(R(rank,rank));
      if (std::abs(R(rank,rank)) < r_min)
	r_min = std::abs(R(rank,rank));
      if (r_min / r_max < rank_threshold)
	break;
    }

    // Extract blocks from Q and R
    R.reshape(rank,rank);
    Q.reshape(Q.numRows(), rank);

    return rank;
  }

  //! Compute column-pivoted QR using modified Gram-Schmidt
  /*!
   * For A an m-by-n matrix, computes A*P = Q*R with R k-by-k upper triangular,
   * Q m-by-k with orthonormal columns, and P an n-by-k permutation matrix.
   * Here k <= min(m,n) is determined by a rank threshold tau provided by the
   * user.  The resulting R will have cond(R) <= 1/tau.  P is returned in the
   * pivot array \c piv and the rank k returned by the function.  Only the first
   * k entries of \c piv will be set.  As with LAPACK, the user can require
   * columns of A to be included in P by setting the corresponding entries
   * of piv to be nonzero on input.  The orthogonality of Q is determined by
   * the weight vector w, defining a weighted inner-product.
   */
  template <typename ordinal_type, typename scalar_type>
  ordinal_type
  CPQR_MGS_threshold(
    const scalar_type& rank_threshold,
    const Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& A,
    const Teuchos::Array<scalar_type>& w,
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& Q,
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& R,
    Teuchos::Array<ordinal_type>& piv) 
  {
    using Teuchos::getCol;
    typedef Teuchos::SerialDenseVector<ordinal_type,scalar_type> SDV;
    typedef Teuchos::SerialDenseMatrix<ordinal_type,scalar_type> SDM;
    ordinal_type m = A.numRows();
    ordinal_type n = A.numCols();
    ordinal_type k = std::min(m,n);
    
    // Make sure Q is the right size
    if (Q.numRows() != m || Q.numCols() != n)
      Q.shape(m,n);
    if (R.numRows() != m || R.numCols() != n)
      R.shape(m,n);
    if (piv.size() != n)
      piv.resize(n);
    Q.assign(A);

    // Compute column norms
    SDV nrm(n);
    for (ordinal_type j=0; j<n; j++) {
      SDV Qj = getCol(Teuchos::View, Q, j);
      nrm(j) = detail::weighted_inner_product(Qj, Qj, w);
    }
    SDV Qtmp(m), Rtmp(m);

    Teuchos::Array<ordinal_type> piv_orig(piv);
    for (ordinal_type i=0; i<n; i++)
      piv[i] = i;

    // Prepivot any columns requested by setting piv[i] != 0
    ordinal_type nfixed = 0;
    for (ordinal_type j=0; j<n; j++) {
      if (piv_orig[j] != 0) {
	if (j != nfixed) {
	  scalar_type tmp = nrm(j);
	  nrm(j) = nrm(nfixed);
	  nrm(nfixed) = tmp;
	  
	  SDV Qpvt = getCol(Teuchos::View, Q, j);
	  SDV Qj = getCol(Teuchos::View, Q, nfixed);
	  Qtmp.assign(Qpvt);
	  Qpvt.assign(Qj);
	  Qj.assign(Qtmp);
	  
	  ordinal_type t = piv[j];
	  piv[j] = piv[nfixed];
	  piv[nfixed] = t;
	}
	++nfixed;
      }
    }
  
    scalar_type sigma = 1.0 + rank_threshold;
    scalar_type r_max = 0.0;
    ordinal_type j=0;
    R.putScalar(0.0);
    while (j < k && sigma >= rank_threshold) {

      SDV Qj = getCol(Teuchos::View, Q, j);

      // Find pivot column if we are passed the pre-pivot stage
      if (j >= nfixed) {
	ordinal_type pvt = j;
	for (ordinal_type l=j+1; l<n; l++)
	  if (nrm(l) > nrm(pvt))
	    pvt = l;

	// Interchange column j and pvt
	if (pvt != j) {
	  SDV Qpvt = getCol(Teuchos::View, Q, pvt);
	  Qtmp.assign(Qpvt);
	  Qpvt.assign(Qj);
	  Qj.assign(Qtmp);

	  SDV Rpvt = getCol(Teuchos::View, R, pvt);
	  SDV Rj = getCol(Teuchos::View, R, j);
	  Rtmp.assign(Rpvt);
	  Rpvt.assign(Rj);
	  Rj.assign(Rtmp);

	  ordinal_type t = piv[pvt];
	  piv[pvt] = piv[j];
	  piv[j] = t;
	}
      }
      
      R(j,j) = std::sqrt(detail::weighted_inner_product(Qj, Qj, w));
      Qj.scale(1.0/R(j,j));
      for (ordinal_type l=j+1; l<n; l++) {
	SDV Ql = getCol(Teuchos::View, Q, l);
	scalar_type t = detail::weighted_inner_product(Qj, Ql, w);
	R(j,l) = t;
	detail::saxpy(1.0, Ql, -t, Qj);

	// Update norms
	nrm(l) = detail::weighted_inner_product(Ql, Ql, w);
      }	

      if (std::abs(R(j,j)) > r_max)
	r_max = R(j,j);
      sigma = std::abs(R(j,j)/r_max);
      j++;
    }

    ordinal_type rank = j;
    if (sigma < rank_threshold)
      rank--;
    Q.reshape(m, rank);
    R.reshape(rank, rank);

    return rank;
  }

  //! Compute column-pivoted QR using modified Gram-Schmidt and reorthogonalization
  /*!
   * For A an m-by-n matrix, computes A*P = Q*R with R k-by-k upper triangular,
   * Q m-by-k with orthonormal columns, and P an n-by-k permutation matrix.
   * Here k <= min(m,n) is determined by a rank threshold tau provided by the
   * user.  The resulting R will have cond(R) <= 1/tau.  P is returned in the
   * pivot array \c piv and the rank k returned by the function.  Only the first
   * k entries of \c piv will be set.  As with LAPACK, the user can require
   * columns of A to be included in P by setting the corresponding entries
   * of piv to be nonzero on input.  The orthogonality of Q is determined by
   * the weight vector w, defining a weighted inner-product.
   */
  template <typename ordinal_type, typename scalar_type>
  ordinal_type
  CPQR_MGS_reorthog_threshold(
    const scalar_type& rank_threshold,
    const Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& A,
    const Teuchos::Array<scalar_type>& w,
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& Q,
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& R,
    Teuchos::Array<ordinal_type>& piv) 
  {
    // First do standard column-pivoted QR
    ordinal_type rank = CPQR_MGS_threshold(rank_threshold, A, w, Q, R, piv);

    // Now apply standard MGS to Q
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type> A2(Q), R2;
    QR_MGS(rank, A2, w, Q, R2);

    return rank;
  }
    
  //! Compute condition number of upper-triangular R
  template <typename ordinal_type, typename scalar_type>
  scalar_type
  cond_R(const Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& R)
  {
    ordinal_type k = R.numRows();
    if (R.numCols() < k)
      k = R.numCols();
    scalar_type r_max = std::abs(R(0,0));
    scalar_type r_min = std::abs(R(0,0));
    for (ordinal_type i=1; i<k; i++) {
      if (std::abs(R(i,i)) > r_max)
	r_max = R(i,i);
      if (std::abs(R(i,i)) < r_min)
	r_min = R(i,i);
    }
    scalar_type cond_r = r_max / r_min;
    return cond_r;
  }
    
  //! Compute weighted QR orthogonalization error
  /*!
   * Computes ||Q^T*W*Q-I||_infinity for Q coming from a weighted QR 
   * factorization.
   */
  template <typename ordinal_type, typename scalar_type>
  scalar_type
  weightedQROrthogonalizationError(
    const Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& Q,
    const Teuchos::Array<scalar_type>& w) 
  {
    ordinal_type m = Q.numRows();
    ordinal_type n = Q.numCols();
    Teuchos::SerialDenseMatrix<ordinal_type, scalar_type> Qt(m,n);
    for (ordinal_type i=0; i<m; i++)
      for (ordinal_type j=0; j<n; j++)
	Qt(i,j) = w[i]*Q(i,j);
    Teuchos::SerialDenseMatrix<ordinal_type, scalar_type> err(n,n);
    err.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Q, Qt, 0.0);
    for (ordinal_type i=0; i<n; i++)
      err(i,i) -= 1.0;
    return err.normInf();
  }

  //! Compute QR orthogonalization error
  /*!
   * Computes ||Q^T*Q-I||_infinity for Q coming from a QR factorization.
   */
  template <typename ordinal_type, typename scalar_type>
  scalar_type
  QROrthogonalizationError(
    const Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& Q) 
  {
    ordinal_type n = Q.numCols();
    Teuchos::SerialDenseMatrix<ordinal_type, scalar_type> err(n,n);
    err.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Q, Q, 0.0);
    for (ordinal_type i=0; i<n; i++)
      err(i,i) -= 1.0;
    return err.normInf();
  }

  //! Compute QR residual error
  /*!
   * Computes ||Q*R-A||_infinity for Q,R coming from QR factorization.
   *
   * Works with thin or full QR, weighted or not.
   */
  template <typename ordinal_type, typename scalar_type>
  scalar_type
  residualQRError(
    const Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>& A,
    const Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>& Q,
    const Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>& R) 
  {
    ordinal_type k = Q.numCols();
    ordinal_type m = Q.numRows();
    Teuchos::SerialDenseMatrix<ordinal_type, scalar_type> AA(
      Teuchos::View, A, m, k, 0, 0);
    Teuchos::SerialDenseMatrix<ordinal_type, scalar_type> err(m,k);
    ordinal_type ret = 
      err.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Q, R, 0.0);
    TEUCHOS_ASSERT(ret == 0);
    err -= AA;
    return err.normInf();
  }

  //! Compute column-pivoted QR residual error
  /*!
   * Computes ||Q*R-A*P||_infinity for Q,R coming from a column-pivoted 
   * QR factorization.
   *
   * Works with thin or full QR, weighted or not.
   */
  template <typename ordinal_type, typename scalar_type>
  scalar_type
  residualCPQRError(
    const Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>& A,
    const Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>& Q,
    const Teuchos::SerialDenseMatrix<ordinal_type, scalar_type>& R,
    const Teuchos::Array<ordinal_type>& piv) 
  {
    ordinal_type k = Q.numCols();
    ordinal_type m = Q.numRows();
    Teuchos::SerialDenseMatrix<ordinal_type, scalar_type> AP(m, k);
    for (ordinal_type j=0; j<k; j++)
      for (ordinal_type i=0; i<m; i++)
	AP(i,j) = A(i,piv[j]);
    Teuchos::SerialDenseMatrix<ordinal_type, scalar_type> err(m,k);
    ordinal_type ret = 
      err.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Q, R, 0.0);
    TEUCHOS_ASSERT(ret == 0);
    err -= AP;

    return err.normInf();
  }

  //! Compute SVD of matrix
  /*!
   * The SVD is computed by the corresponding LAPACK function.
   */
  template <typename ordinal_type, typename scalar_type>
  void
  svd(const Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& A,
      Teuchos::Array<scalar_type>& s,
      Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& U,
      Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& Vt) 
  {
    Teuchos::LAPACK<ordinal_type,scalar_type> lapack;
    ordinal_type m = A.numRows();
    ordinal_type n = A.numCols();
    ordinal_type k = std::min(m,n);
    ordinal_type lda = A.stride();

    // Copy A since GESVD overwrites it (always)
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type> AA(
      Teuchos::Copy, A, m, n);

    // Resize appropriately
    if (U.numRows() != m || U.numCols() != m)
      U.shape(m,m);
    if (Vt.numRows() != n || Vt.numCols() != n)
      Vt.shape(n,n);
    if (s.size() != k)
      s.resize(k);
    ordinal_type ldu = U.stride();
    ordinal_type ldv = Vt.stride();
    
    // Workspace query
    Teuchos::Array<scalar_type> work(1);
    ordinal_type lwork = -1;
    ordinal_type info;
    lapack.GESVD('A', 'A', m, n, AA.values(), lda, &s[0], U.values(), ldu, 
		 Vt.values(), ldv, &work[0], lwork, NULL, &info);
    TEUCHOS_TEST_FOR_EXCEPTION(
      info < 0, std::logic_error, "dgesvd returned info = " << info);

    // Do SVD
    lwork = work[0];
    work.resize(lwork);
    lapack.GESVD('A', 'A', m, n, AA.values(), lda, &s[0], U.values(), ldu, 
		 Vt.values(), ldv, &work[0], lwork, NULL, &info);
    TEUCHOS_TEST_FOR_EXCEPTION(
      info < 0, std::logic_error, "dgesvd returned info = " << info);
    
  }

  template <typename ordinal_type, typename scalar_type>
  ordinal_type
  svd_threshold(
    const scalar_type& rank_threshold,
    const Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& A,
    Teuchos::Array<scalar_type>& s,
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& U,
    Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& Vt) 
  {
    // Compute full SVD
    svd(A, s, U, Vt);

    // Find leading block of S such that cond(S) <= 1/tau
    ordinal_type rank = 0;
    ordinal_type m = s.size();
    while (rank < m && s[rank]/s[0] > rank_threshold) rank++;

    // Extract blocks from s, U and Vt
    s.resize(rank);
    U.reshape(U.numRows(),rank);
    Vt.reshape(rank, Vt.numCols());

    return rank;
  }

}

#endif
