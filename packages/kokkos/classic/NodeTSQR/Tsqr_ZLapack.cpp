//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include <Tsqr_Lapack.hpp>
#include <complex>


extern "C" void F77_BLAS_MANGLE(zlarnv, ZLARNV)
  (const int* const IDIST,
   int ISEED[],
   const int* const N,
   std::complex<double> X[]);

extern "C" void F77_BLAS_MANGLE(zpotri, ZPOTRI)
  (const char* const UPLO,
   const int* const N,
   std::complex<double> A[],
   const int* const LDA,
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(zpotrf, ZPOTRF)
  (const char* const UPLO,
   const int* const N,
   std::complex<double> A[],
   const int* const LDA,
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(zpotrs, ZPOTRS)
  (const char* const UPLO,
   const int* const N,
   const int* const NRHS,
   const std::complex<double> A[],
   const int* const LDA,
   std::complex<double> B[],
   const int* const LDB,
   int* const INFO);

#ifdef HAVE_LAPACK_ZLARFGP
extern "C" void F77_BLAS_MANGLE(zlarfgp,ZLARFGP)
  (const int* const N,    // IN
   std::complex<double>* const ALPHA,   // IN/OUT
   std::complex<double> X[],            // IN/OUT
   const int* const INCX, // IN
   std::complex<double>* const TAU);    // OUT
#else
#  ifdef HAVE_LAPACK_ZLARFP
extern "C" void F77_BLAS_MANGLE(zlarfp,ZLARFP)
  (const int* const N,    // IN
   std::complex<double>* const ALPHA,   // IN/OUT
   std::complex<double> X[],            // IN/OUT
   const int* const INCX, // IN
   std::complex<double>* const TAU);    // OUT
#  else
extern "C" void F77_BLAS_MANGLE(zlarfg,ZLARFG)
  (const int* const N,    // IN
   std::complex<double>* const ALPHA,   // IN/OUT
   std::complex<double> X[],            // IN/OUT
   const int* const INCX, // IN
   std::complex<double>* const TAU);    // OUT
#  endif // HAVE_LAPACK_ZLARFP
#endif // HAVE_LAPACK_ZLARFGP

extern "C" void F77_BLAS_MANGLE(zgeqrf, ZGEQRF)
  (const int* const M,
   const int* const N,
   std::complex<double> A[],
   const int* const LDA,
   std::complex<double> TAU[],
   std::complex<double> WORK[],
   const int* const LWORK,
   int* const INFO);

#ifdef HAVE_LAPACK_ZGEQRFP
extern "C" void F77_BLAS_MANGLE(zgeqrfp, ZGEQRFP)
  (const int* const M,
   const int* const N,
   std::complex<double> A[],
   const int* const LDA,
   std::complex<double> TAU[],
   std::complex<double> WORK[],
   const int* const LWORK,
   int* const INFO);
#endif // HAVE_LAPACK_ZGEQRFP

extern "C" void F77_BLAS_MANGLE(zgeqr2, ZGEQR2)
  (const int* const M,
   const int* const N,
   std::complex<double> A[],
   const int* const LDA,
   std::complex<double> TAU[],
   std::complex<double> WORK[],
   int* const INFO);

#ifdef HAVE_LAPACK_ZGEQR2P
extern "C" void F77_BLAS_MANGLE(zgeqr2p, ZGEQR2P)
  (const int* const M,
   const int* const N,
   std::complex<double> A[],
   const int* const LDA,
   std::complex<double> TAU[],
   std::complex<double> WORK[],
   int* const INFO);
#endif // HAVE_LAPACK_ZGEQR2P

// In the complex case, Q is called UNitary rather than ORthogonal.
// This is why we have ZUNGQR and CUNGQR, rather than ZORGQR and
// CORGQR.  The interface is exactly the same as in the real case,
// though, so our LAPACK::ORMQR(), etc. wrappers have the same name
// for both the real and the complex cases.

extern "C" void F77_BLAS_MANGLE(zungqr, ZUNGQR)
  (const int* const M,
   const int* const N,
   const int* const K,
   std::complex<double> A[],
   const int* const LDA,
   std::complex<double> TAU[],
   std::complex<double> WORK[],
   const int* const LWORK,
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(zunmqr, ZUNMQR)
  (const char* const SIDE,
   const char* const TRANS,
   const int* const M,
   const int* const N,
   const int* const K,
   const std::complex<double> A[],
   const int* const LDA,
   const std::complex<double> TAU[],
   std::complex<double> C[],
   const int* const LDC,
   std::complex<double> WORK[],
   const int* const LWORK,
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(zunm2r, ZUNM2R)
  (const char* const SIDE,
   const char* const TRANS,
   const int* const M,
   const int* const N,
   const int* const K,
   const std::complex<double> A[],
   const int* const LDA,
   const std::complex<double> TAU[],
   std::complex<double> C[],
   const int* const LDC,
   std::complex<double> WORK[],
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(zgesvd, ZGESVD) 
  (const char* const JOBU, 
   const char* const JOBVT, 
   const int* const M, 
   const int* const N, 
   std::complex<double> A[], 
   const int* const LDA,
   double S[], 
   std::complex<double> U[], 
   const int* const LDU, 
   std::complex<double> VT[], 
   const int* const LDVT, 
   std::complex<double> work[],
   const int* const LWORK,
   double RWORK[],
   int* const INFO);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  // If _GEQRFP is available, LAPACK::GEQRF() calls it.  If _LARFP is
  // available, LAPACK::GEQRF() calls _GEQRF, which uses _LARFP.
#ifdef HAVE_LAPACK_ZGEQRFP
  template <>
  bool LAPACK<int, std::complex<double> >::QR_produces_R_factor_with_nonnegative_diagonal() { return true; }
#else
#  ifdef HAVE_LAPACK_ZLARFP
  template <>
  bool LAPACK<int, std::complex<double> >::QR_produces_R_factor_with_nonnegative_diagonal() { return true; }
#  else
  template <>
  bool LAPACK<int, std::complex<double> >::QR_produces_R_factor_with_nonnegative_diagonal() { return false; }
#  endif
#endif

  ////////////////////////////////////////////////////////////////////////////
  // LARFP (implemented with _LARFGP if available, else with _LARFP if
  // available, else fall back to _LARFG)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void 
  LAPACK<int, std::complex<double> >::
  LARFP (const int n, 
	 std::complex<double>& alpha, 
	 std::complex<double> x[], 
	 const int incx, 
	 std::complex<double>& tau)
  {
#ifdef HAVE_LAPACK_ZLARFGP
    F77_BLAS_MANGLE(zlarfgp,ZLARFGP) (&n, &alpha, x, &incx, &tau);
#else // Don't HAVE_LAPACK_CLARFGP
#  ifdef HAVE_LAPACK_ZLARFP
    F77_BLAS_MANGLE(zlarfp,ZLARFP) (&n, &alpha, x, &incx, &tau);
#  else
    F77_BLAS_MANGLE(zlarfg,ZLARFG) (&n, &alpha, x, &incx, &tau);
#  endif // HAVE_LAPACK_ZLARFP
#endif // HAVE_LAPACK_ZLARFGP
  }

  ////////////////////////////////////////////////////////////////////////////
  // GEQRF (implemented with _GEQRFP if available, else fall back to _GEQRF)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void
  LAPACK<int, std::complex<double> >::
  GEQRF (const int m,
	 const int n, 
	 std::complex<double> A[],
	 const int lda, 
	 std::complex<double> tau[],
	 std::complex<double> work[],
	 const int lwork,
	 int* const INFO)
  {
#ifdef HAVE_LAPACK_ZGEQRFP
    F77_BLAS_MANGLE(zgeqrfp, ZGEQRFP) (&m, &n, A, &lda, tau, work, &lwork, INFO);
#else
    F77_BLAS_MANGLE(zgeqrf, ZGEQRF) (&m, &n, A, &lda, tau, work, &lwork, INFO);
#endif // HAVE_LAPACK_ZGEQRFP
  }

  ////////////////////////////////////////////////////////////////////////////
  // GEQR2 (implemented with _GEQR2P if available, else fall back to _GEQR2)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void
  LAPACK<int, std::complex<double> >::
  GEQR2 (const int m,
	 const int n, 
	 std::complex<double> A[],
	 const int lda, 
	 std::complex<double> tau[],
	 std::complex<double> work[],
	 int* const INFO)
  {
#ifdef HAVE_LAPACK_ZGEQR2P
    F77_BLAS_MANGLE(zgeqr2p, ZGEQR2P) (&m, &n, A, &lda, tau, work, INFO);
#else
    F77_BLAS_MANGLE(zgeqr2, ZGEQR2) (&m, &n, A, &lda, tau, work, INFO);
#endif // HAVE_LAPACK_ZGEQR2P
  }

  template <>
  void
  LAPACK<int, std::complex<double> >::
  ORMQR (const char* const side,
	 const char* const trans,
	 const int m,
	 const int n,
	 const int k,
	 const std::complex<double> A[],
	 const int lda,
	 const std::complex<double> tau[],
	 std::complex<double> C[],
	 const int ldc,
	 std::complex<double> work[],
	 const int lwork,
	 int* const INFO)
  {
    F77_BLAS_MANGLE(zunmqr, ZUNMQR) 
      (side, trans, &m, &n, &k, A, &lda, tau, C, &ldc, work, &lwork, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<double> >::
  ORM2R (const char* const side,
	 const char* const trans,
	 const int m,
	 const int n,
	 const int k,
	 const std::complex<double> A[],
	 const int lda,
	 const std::complex<double> tau[],
	 std::complex<double> C[],
	 const int ldc,
	 std::complex<double> work[],
	 int* const INFO)
  {
    F77_BLAS_MANGLE(zunm2r, ZUNM2R) 
      (side, trans, &m, &n, &k, A, &lda, tau, C, &ldc, work, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<double> >::
  ORGQR (const int m,
	 const int n,
	 const int k,
	 std::complex<double> A[],
	 const int lda,
	 std::complex<double> tau[],
	 std::complex<double> work[],
	 const int lwork,
	 int* const INFO)
  {
    F77_BLAS_MANGLE(zungqr, ZUNGQR) 
      (&m, &n, &k, A, &lda, tau, work, &lwork, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<double> >::
  POTRF (const char* const uplo,
	 const int n,
	 std::complex<double> A[],
	 const int lda,
	 int* const INFO)
  {
    F77_BLAS_MANGLE(zpotrf, ZPOTRF) (uplo, &n, A, &lda, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<double> >::
  POTRS (const char* const uplo,
	 const int n,
	 const int nrhs,
	 const std::complex<double> A[],
	 const int lda,
	 std::complex<double> B[],
	 const int ldb,
	 int* const INFO)
  {
    F77_BLAS_MANGLE(zpotrs, ZPOTRS) 
      (uplo, &n, &nrhs, A, &lda, B, &ldb, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<double> >::
  POTRI (const char* const uplo, 
	 const int n, 
	 std::complex<double> A[], 
	 const int lda, 
	 int* const INFO)
  {
    F77_BLAS_MANGLE(zpotri, ZPOTRI) (uplo, &n, A, &lda, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<double> >::
  LARNV (const int idist, 
	 int iseed[],
	 const int n,
	 std::complex<double> x[])
  {
    F77_BLAS_MANGLE(zlarnv, ZLARNV) (&idist, iseed, &n, x);
  }

  template <>
  void
  LAPACK<int, std::complex<double> >::
  GESVD (const char* const jobu,
	 const char* const jobvt,
	 const int m,
	 const int n,
	 std::complex<double> A[],
	 const int lda,
	 double s[],
	 std::complex<double> U[],
	 const int ldu,
	 std::complex<double> VT[],
	 const int ldvt,
	 std::complex<double> work[],
	 const int lwork,
	 double rwork[],
	 int* const INFO)
  {
    F77_BLAS_MANGLE(zgesvd, ZGESVD) (jobu, jobvt, &m, &n, 
				     A, &lda, s, 
				     U, &ldu, VT, &ldvt, 
				     work, &lwork, rwork, INFO);
  }


} // namespace TSQR
