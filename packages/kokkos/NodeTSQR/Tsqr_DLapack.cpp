//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
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
// ************************************************************************
//@HEADER

#include <Tsqr_Lapack.hpp>

extern "C" void F77_BLAS_MANGLE(dlarnv, DLARNV)
  (const int* const IDIST,
   int ISEED[],
   const int* const N,
   double X[]);

extern "C" void F77_BLAS_MANGLE(dpotri, DPOTRI)
  (const char* const UPLO,
   const int* const N,
   double A[],
   const int* const LDA,
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(dpotrf, DPOTRF)
  (const char* const UPLO,
   const int* const N,
   double A[],
   const int* const LDA,
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(dpotrs, DPOTRS)
  (const char* const UPLO,
   const int* const N,
   const int* const NRHS,
   const double A[],
   const int* const LDA,
   double B[],
   const int* const LDB,
   int* const INFO);

#ifdef HAVE_LAPACK_DLARFGP
extern "C" void F77_BLAS_MANGLE(dlarfgp,DLARFGP)
  (const int* const N,    // IN
   double* const ALPHA,   // IN/OUT
   double X[],            // IN/OUT
   const int* const INCX, // IN
   double* const TAU);    // OUT
#else
#  ifdef HAVE_LAPACK_DLARFP
extern "C" void F77_BLAS_MANGLE(dlarfp,DLARFP)
  (const int* const N,    // IN
   double* const ALPHA,   // IN/OUT
   double X[],            // IN/OUT
   const int* const INCX, // IN
   double* const TAU);    // OUT
#  else
extern "C" void F77_BLAS_MANGLE(dlarfg,DLARFG)
  (const int* const N,    // IN
   double* const ALPHA,   // IN/OUT
   double X[],            // IN/OUT
   const int* const INCX, // IN
   double* const TAU);    // OUT
#  endif // HAVE_LAPACK_DLARFP
#endif // HAVE_LAPACK_DLARFGP

extern "C" void F77_BLAS_MANGLE(dgeqrf, DGEQRF)
  (const int* const M,
   const int* const N,
   double A[],
   const int* const LDA,
   double TAU[],
   double WORK[],
   const int* const LWORK,
   int* const INFO);

#ifdef HAVE_LAPACK_DGEQRFP
extern "C" void F77_BLAS_MANGLE(dgeqrfp, DGEQRFP)
  (const int* const M,
   const int* const N,
   double A[],
   const int* const LDA,
   double TAU[],
   double WORK[],
   const int* const LWORK,
   int* const INFO);
#endif // HAVE_LAPACK_DGEQRFP

extern "C" void F77_BLAS_MANGLE(dgeqr2, DGEQR2)
  (const int* const M,
   const int* const N,
   double A[],
   const int* const LDA,
   double TAU[],
   double WORK[],
   int* const INFO);

#ifdef HAVE_LAPACK_DGEQR2P
extern "C" void F77_BLAS_MANGLE(dgeqr2p, DGEQR2P)
  (const int* const M,
   const int* const N,
   double A[],
   const int* const LDA,
   double TAU[],
   double WORK[],
   int* const INFO);
#endif // HAVE_LAPACK_DGEQR2P

extern "C" void F77_BLAS_MANGLE(dormqr, DORMQR)
  (const char* const SIDE,
   const char* const TRANS,
   const int* const M,
   const int* const N,
   const int* const K,
   const double A[],
   const int* const LDA,
   const double TAU[],
   double C[],
   const int* const LDC,
   double WORK[],
   const int* const LWORK,
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(dorm2r, DORM2R)
  (const char* const SIDE,
   const char* const TRANS,
   const int* const M,
   const int* const N,
   const int* const K,
   const double A[],
   const int* const LDA,
   const double TAU[],
   double C[],
   const int* const LDC,
   double WORK[],
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(dorgqr, DORGQR)
  (const int* const M,
   const int* const N,
   const int* const K,
   double A[],
   const int* const LDA,
   double TAU[],
   double WORK[],
   const int* const LWORK,
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(dgesvd, DGESVD) 
  (const char* const JOBU, 
   const char* const JOBVT, 
   const int* const M, 
   const int* const N, 
   double A[], 
   const int* const LDA,
   double S[], 
   double U[], 
   const int* const LDU, 
   double VT[], 
   const int* const LDVT, 
   double work[],
   const int* const LWORK,
   double RWORK[],
   int* const INFO);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  // If _GEQRFP is available, LAPACK::GEQRF() calls it.  If _LARFP is
  // available, LAPACK::GEQRF() calls _GEQRF, which uses _LARFP.
#ifdef HAVE_LAPACK_DGEQRFP
  template <>
  bool LAPACK<int, double >::QR_produces_R_factor_with_nonnegative_diagonal() { return true; }
#else
#  ifdef HAVE_LAPACK_DLARFP
  template <>
  bool LAPACK<int, double >::QR_produces_R_factor_with_nonnegative_diagonal() { return true; }
#  else
  template <>
  bool LAPACK<int, double >::QR_produces_R_factor_with_nonnegative_diagonal() { return false; }
#  endif
#endif

  ////////////////////////////////////////////////////////////////////////////
  // LARFP (implemented with _LARFGP if available, else with _LARFP if
  // available, else fall back to _LARFG)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void 
  LAPACK<int, double >::LARFP (const int n, 
			       double& alpha, 
			       double x[], 
			       const int incx, 
			       double& tau)
  {
#ifdef HAVE_LAPACK_DLARFGP
    F77_BLAS_MANGLE(dlarfgp,DLARFGP) (&n, &alpha, x, &incx, &tau);
#else // Don't HAVE_LAPACK_DLARFGP
#  ifdef HAVE_LAPACK_DLARFP
    F77_BLAS_MANGLE(dlarfp,DLARFP) (&n, &alpha, x, &incx, &tau);
#  else
    F77_BLAS_MANGLE(dlarfg,DLARFG) (&n, &alpha, x, &incx, &tau);
#  endif // HAVE_LAPACK_DLARFP
#endif // HAVE_LAPACK_DLARFGP
  }

  ////////////////////////////////////////////////////////////////////////////
  // GEQRF (implemented with _GEQRFP if available, else fall back to _GEQRF)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void
  LAPACK<int, double >::GEQRF (const int m,
			       const int n, 
			       double A[],
			       const int lda, 
			       double tau[],
			       double work[],
			       const int lwork,
			       int* const INFO)
  {
#ifdef HAVE_LAPACK_DGEQRFP
    F77_BLAS_MANGLE(dgeqrfp, DGEQRFP) 
      (&m, &n, A, &lda, tau, work, &lwork, INFO);
#else
    F77_BLAS_MANGLE(dgeqrf, DGEQRF) 
      (&m, &n, A, &lda, tau, work, &lwork, INFO);
#endif // HAVE_LAPACK_DGEQRFP
  }

  ////////////////////////////////////////////////////////////////////////////
  // GEQR2 (implemented with _GEQR2P if available, else fall back to _GEQR2)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void
  LAPACK<int, double >::GEQR2 (const int m,
			       const int n, 
			       double A[],
			       const int lda, 
			       double tau[],
			       double work[],
			       int* const INFO)
  {
#ifdef HAVE_LAPACK_DGEQR2P
    F77_BLAS_MANGLE(dgeqr2p, DGEQR2P) (&m, &n, A, &lda, tau, work, INFO);
#else
    F77_BLAS_MANGLE(dgeqr2, DGEQR2) (&m, &n, A, &lda, tau, work, INFO);
#endif // HAVE_LAPACK_DGEQR2P
  }

  template <>
  void
  LAPACK<int, double >::ORMQR (const char* const side,
			       const char* const trans,
			       const int m,
			       const int n,
			       const int k,
			       const double A[],
			       const int lda,
			       const double tau[],
			       double C[],
			       const int ldc,
			       double work[],
			       const int lwork,
			       int* const INFO)
  {
    F77_BLAS_MANGLE(dormqr, DORMQR) 
      (side, trans, &m, &n, &k, A, &lda, tau, C, &ldc, work, &lwork, INFO);
  }

  template <>
  void
  LAPACK<int, double >::ORM2R (const char* const side,
			       const char* const trans,
			       const int m,
			       const int n,
			       const int k,
			       const double A[],
			       const int lda,
			       const double tau[],
			       double C[],
			       const int ldc,
			       double work[],
			       int* const INFO)
  {
    F77_BLAS_MANGLE(dorm2r, DORM2R) 
      (side, trans, &m, &n, &k, A, &lda, tau, C, &ldc, work, INFO);
  }

  template <>
  void
  LAPACK<int, double >::ORGQR (const int m,
			       const int n,
			       const int k,
			       double A[],
			       const int lda,
			       double tau[],
			       double work[],
			       const int lwork,
			       int* const INFO)
  {
    F77_BLAS_MANGLE(dorgqr, DORGQR) 
      (&m, &n, &k, A, &lda, tau, work, &lwork, INFO);
  }

  template <>
  void
  LAPACK<int, double >::POTRF (const char* const uplo,
			       const int n,
			       double A[],
			       const int lda,
			       int* const INFO)
  {
    F77_BLAS_MANGLE(dpotrf, DPOTRF) (uplo, &n, A, &lda, INFO);
  }

  template <>
  void
  LAPACK<int, double >::POTRS (const char* const uplo,
			       const int n,
			       const int nrhs,
			       const double A[],
			       const int lda,
			       double B[],
			       const int ldb,
			       int* const INFO)
  {
    F77_BLAS_MANGLE(dpotrs, DPOTRS) 
      (uplo, &n, &nrhs, A, &lda, B, &ldb, INFO);
  }

  template <>
  void
  LAPACK<int, double >::POTRI (const char* const uplo, 
			       const int n, 
			       double A[], 
			       const int lda, 
			       int* const INFO)
  {
    F77_BLAS_MANGLE(dpotri, DPOTRI) (uplo, &n, A, &lda, INFO);
  }

  template <>
  void
  LAPACK<int, double >::LARNV (const int idist, 
			       int iseed[],
			       const int n,
			       double x[])
  {
    F77_BLAS_MANGLE(dlarnv, DLARNV) (&idist, iseed, &n, x);
  }

  template <>
  void
  LAPACK<int, double >::GESVD (const char* const jobu,
			       const char* const jobvt,
			       const int m,
			       const int n,
			       double A[],
			       const int lda,
			       double s[],
			       double U[],
			       const int ldu,
			       double VT[],
			       const int ldvt,
			       double work[],
			       const int lwork,
			       double rwork[],
			       int* const INFO)
  {
    F77_BLAS_MANGLE(dgesvd, DGESVD) (jobu, jobvt, &m, &n, 
				     A, &lda, s, 
				     U, &ldu, VT, &ldvt, 
				     work, &lwork, rwork, INFO);
  }

} // namespace TSQR
