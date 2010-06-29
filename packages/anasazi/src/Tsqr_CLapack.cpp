#include <Tsqr_Lapack.hpp>
#include <Tsqr_FortranCInterface.hpp>
#include <Tsqr_Config.hpp>
#include <complex>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

extern "C" void FortranCInterface_GLOBAL(cpotri, CPOTRI)
  (const char* const UPLO,
   const int* const N,
   std::complex<float> A[],
   const int* const LDA,
   int* const INFO);

extern "C" void FortranCInterface_GLOBAL(cpotrf, CPOTRF)
  (const char* const UPLO,
   const int* const N,
   std::complex<float> A[],
   const int* const LDA,
   int* const INFO);

extern "C" void FortranCInterface_GLOBAL(cpotrs, CPOTRS)
  (const char* const UPLO,
   const int* const N,
   const int* const NRHS,
   const std::complex<float> A[],
   const int* const LDA,
   std::complex<float> B[],
   const int* const LDB,
   int* const INFO);

#ifdef HAVE_LAPACK_CLARFGP
extern "C" void FortranCInterface_GLOBAL(clarfgp,CLARFGP)
  (const int* const N,                 // IN
   std::complex<float>* const ALPHA,   // IN/OUT
   std::complex<float> X[],            // IN/OUT
   const int* const INCX,              // IN
   std::complex<float>* const TAU);    // OUT
#else
#  ifdef HAVE_LAPACK_CLARFP
extern "C" void FortranCInterface_GLOBAL(clarfp,CLARFP)
  (const int* const N,                 // IN
   std::complex<float>* const ALPHA,   // IN/OUT
   std::complex<float> X[],            // IN/OUT
   const int* const INCX,              // IN
   std::complex<float>* const TAU);    // OUT
#  else
extern "C" void FortranCInterface_GLOBAL(clarfg,CLARFG)
  (const int* const N,                 // IN
   std::complex<float>* const ALPHA,   // IN/OUT
   std::complex<float> X[],            // IN/OUT
   const int* const INCX,              // IN
   std::complex<float>* const TAU);    // OUT
#  endif // HAVE_LAPACK_CLARFP
#endif // HAVE_LAPACK_CLARFGP

extern "C" void FortranCInterface_GLOBAL(cgeqrf, CGEQRF)
  (const int* const M,
   const int* const N,
   std::complex<float> A[],
   const int* const LDA,
   std::complex<float> TAU[],
   std::complex<float> WORK[],
   const int* const LWORK,
   int* const INFO);

#ifdef HAVE_LAPACK_CGEQRFP
extern "C" void FortranCInterface_GLOBAL(cgeqrfp, CGEQRFP)
  (const int* const M,
   const int* const N,
   std::complex<float> A[],
   const int* const LDA,
   std::complex<float> TAU[],
   std::complex<float> WORK[],
   const int* const LWORK,
   int* const INFO);
#endif // HAVE_LAPACK_CGEQRFP

extern "C" void FortranCInterface_GLOBAL(cgeqr2, CGEQR2)
  (const int* const M,
   const int* const N,
   std::complex<float> A[],
   const int* const LDA,
   std::complex<float> TAU[],
   std::complex<float> WORK[],
   int* const INFO);

#ifdef HAVE_LAPACK_CGEQR2P
extern "C" void FortranCInterface_GLOBAL(cgeqr2p, CGEQR2P)
  (const int* const M,
   const int* const N,
   std::complex<float> A[],
   const int* const LDA,
   std::complex<float> TAU[],
   std::complex<float> WORK[],
   int* const INFO);
#endif // HAVE_LAPACK_CGEQR2P

// In the complex case, Q is called UNitary rather than ORthogonal.
// This is why we have ZUNGQR and CUNGQR, rather than ZORGQR and
// CORGQR.  The interface is exactly the same as in the real case,
// though, so our LAPACK::ORMQR(), etc. wrappers have the same name
// for both the real and the complex cases.

extern "C" void FortranCInterface_GLOBAL(cungqr, CUNGQR)
  (const int* const M,
   const int* const N,
   const int* const K,
   std::complex<float> A[],
   const int* const LDA,
   std::complex<float> TAU[],
   std::complex<float> WORK[],
   const int* const LWORK,
   int* const INFO);

extern "C" void FortranCInterface_GLOBAL(cunm2r, CUNM2R)
  (const char* const SIDE,
   const char* const TRANS,
   const int* const M,
   const int* const N,
   const int* const K,
   const std::complex<float> A[],
   const int* const LDA,
   const std::complex<float> TAU[],
   std::complex<float> C[],
   const int* const LDC,
   std::complex<float> WORK[],
   int* const INFO);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  // If _GEQRFP is available, LAPACK::GEQRF() calls it.  If _LARFP is
  // available, LAPACK::GEQRF() calls _GEQRF, which uses _LARFP.
#ifdef HAVE_LAPACK_CGEQRFP
  template <>
  const bool LAPACK<int, std::complex<float> >::QR_produces_R_factor_with_nonnegative_diagonal = true;
#else
#  ifdef HAVE_LAPACK_CLARFP
  template <>
  const bool LAPACK<int, std::complex<float> >::QR_produces_R_factor_with_nonnegative_diagonal = true;
#  else
  template <>
  const bool LAPACK<int, std::complex<float> >::QR_produces_R_factor_with_nonnegative_diagonal = false;
#  endif
#endif

  ////////////////////////////////////////////////////////////////////////////
  // LARFP (implemented with _LARFGP if available, else with _LARFP if
  // available, else fall back to _LARFG)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void 
  LAPACK<int, std::complex<float> >::LARFP (const int n, 
					    std::complex<float>& alpha, 
					    std::complex<float> x[], 
					    const int incx, 
					    std::complex<float>& tau)
  {
#ifdef HAVE_LAPACK_CLARFGP
    FortranCInterface_GLOBAL(clarfgp,CLARFGP) (&n, &alpha, x, &incx, &tau);
#  ifdef HAVE_LAPACK_CLARFP
    FortranCInterface_GLOBAL(clarfp,CLARFP) (&n, &alpha, x, &incx, &tau);
#  else
    FortranCInterface_GLOBAL(clarfg,CLARFG) (&n, &alpha, x, &incx, &tau);
#  endif // HAVE_LAPACK_CLARFP
#endif // HAVE_LAPACK_CLARFGP
  }

  ////////////////////////////////////////////////////////////////////////////
  // GEQRF (implemented with _GEQRFP if available, else fall back to _GEQRF)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void
  LAPACK<int, std::complex<float> >::GEQRF (const int m,
					    const int n, 
					    std::complex<float> A[],
					    const int lda, 
					    std::complex<float> tau[],
					    std::complex<float> work[],
					    const int lwork,
					    int* const INFO)
  {
#ifdef HAVE_LAPACK_CGEQRFP
    FortranCInterface_GLOBAL(cgeqrfp, CGEQRFP) (&m, &n, A, &lda, tau, work, &lwork, INFO);
#else
    FortranCInterface_GLOBAL(cgeqrf, CGEQRF) (&m, &n, A, &lda, tau, work, &lwork, INFO);
#endif // HAVE_LAPACK_CGEQRFP
  }

  ////////////////////////////////////////////////////////////////////////////
  // GEQR2 (implemented with _GEQR2P if available, else fall back to _GEQR2)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void
  LAPACK<int, std::complex<float> >::GEQR2 (const int m,
					    const int n, 
					    std::complex<float> A[],
					    const int lda, 
					    std::complex<float> tau[],
					    std::complex<float> work[],
					    int* const INFO)
  {
#ifdef HAVE_LAPACK_CGEQR2P
    FortranCInterface_GLOBAL(cgeqr2p, CGEQR2P) (&m, &n, A, &lda, tau, work, INFO);
#else
    FortranCInterface_GLOBAL(cgeqr2, CGEQR2) (&m, &n, A, &lda, tau, work, INFO);
#endif // HAVE_LAPACK_CGEQR2P
  }

  template <>
  void
  LAPACK<int, std::complex<float> >::ORM2R (const char* const side,
					    const char* const trans,
					    const int m,
					    const int n,
					    const int k,
					    const std::complex<float> A[],
					    const int lda,
					    const std::complex<float> tau[],
					    std::complex<float> C[],
					    const int ldc,
					    std::complex<float> work[],
					    int* const INFO)
  {
    FortranCInterface_GLOBAL(cunm2r, CUNM2R) (side, trans, &m, &n, &k, A, &lda, tau, C, &ldc, work, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<float> >::ORGQR (const int m,
					    const int n,
					    const int k,
					    std::complex<float> A[],
					    const int lda,
					    std::complex<float> tau[],
					    std::complex<float> work[],
					    const int lwork,
					    int* const INFO)
  {
    FortranCInterface_GLOBAL(cungqr, CUNGQR) (&m, &n, &k, A, &lda, tau, work, &lwork, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<float> >::POTRF (const char* const uplo,
					    const int n,
					    std::complex<float> A[],
					    const int lda,
					    int* const INFO)
  {
    FortranCInterface_GLOBAL(cpotrf, CPOTRF) (uplo, &n, A, &lda, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<float> >::POTRS (const char* const uplo,
					    const int n,
					    const int nrhs,
					    const std::complex<float> A[],
					    const int lda,
					    std::complex<float> B[],
					    const int ldb,
					    int* const INFO)
  {
    FortranCInterface_GLOBAL(cpotrs, CPOTRS) (uplo, &n, &nrhs, A, &lda, B, &ldb, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<float> >::POTRI (const char* const uplo, 
					    const int n, 
					    std::complex<float> A[], 
					    const int lda, 
					    int* const INFO)
  {
    FortranCInterface_GLOBAL(cpotri, CPOTRI) (uplo, &n, A, &lda, INFO);
  }

} // namespace TSQR
