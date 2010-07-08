#include <Tsqr_Lapack.hpp>
#include <Tsqr_FortranCInterface.hpp>
#include <Tsqr_Config.hpp>
#include <complex>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

extern "C" void FortranCInterface_GLOBAL(zlarnv, ZLARNV)
  (const int* const IDIST,
   int ISEED[],
   const int* const N,
   std::complex<double> X[]);

extern "C" void FortranCInterface_GLOBAL(zpotri, ZPOTRI)
  (const char* const UPLO,
   const int* const N,
   std::complex<double> A[],
   const int* const LDA,
   int* const INFO);

extern "C" void FortranCInterface_GLOBAL(zpotrf, ZPOTRF)
  (const char* const UPLO,
   const int* const N,
   std::complex<double> A[],
   const int* const LDA,
   int* const INFO);

extern "C" void FortranCInterface_GLOBAL(zpotrs, ZPOTRS)
  (const char* const UPLO,
   const int* const N,
   const int* const NRHS,
   const std::complex<double> A[],
   const int* const LDA,
   std::complex<double> B[],
   const int* const LDB,
   int* const INFO);

#ifdef HAVE_LAPACK_ZLARFGP
extern "C" void FortranCInterface_GLOBAL(zlarfgp,ZLARFGP)
  (const int* const N,    // IN
   std::complex<double>* const ALPHA,   // IN/OUT
   std::complex<double> X[],            // IN/OUT
   const int* const INCX, // IN
   std::complex<double>* const TAU);    // OUT
#else
#  ifdef HAVE_LAPACK_ZLARFP
extern "C" void FortranCInterface_GLOBAL(zlarfp,ZLARFP)
  (const int* const N,    // IN
   std::complex<double>* const ALPHA,   // IN/OUT
   std::complex<double> X[],            // IN/OUT
   const int* const INCX, // IN
   std::complex<double>* const TAU);    // OUT
#  else
extern "C" void FortranCInterface_GLOBAL(zlarfg,ZLARFG)
  (const int* const N,    // IN
   std::complex<double>* const ALPHA,   // IN/OUT
   std::complex<double> X[],            // IN/OUT
   const int* const INCX, // IN
   std::complex<double>* const TAU);    // OUT
#  endif // HAVE_LAPACK_ZLARFP
#endif // HAVE_LAPACK_ZLARFGP

extern "C" void FortranCInterface_GLOBAL(zgeqrf, ZGEQRF)
  (const int* const M,
   const int* const N,
   std::complex<double> A[],
   const int* const LDA,
   std::complex<double> TAU[],
   std::complex<double> WORK[],
   const int* const LWORK,
   int* const INFO);

#ifdef HAVE_LAPACK_ZGEQRFP
extern "C" void FortranCInterface_GLOBAL(zgeqrfp, ZGEQRFP)
  (const int* const M,
   const int* const N,
   std::complex<double> A[],
   const int* const LDA,
   std::complex<double> TAU[],
   std::complex<double> WORK[],
   const int* const LWORK,
   int* const INFO);
#endif // HAVE_LAPACK_ZGEQRFP

extern "C" void FortranCInterface_GLOBAL(zgeqr2, ZGEQR2)
  (const int* const M,
   const int* const N,
   std::complex<double> A[],
   const int* const LDA,
   std::complex<double> TAU[],
   std::complex<double> WORK[],
   int* const INFO);

#ifdef HAVE_LAPACK_ZGEQR2P
extern "C" void FortranCInterface_GLOBAL(zgeqr2p, ZGEQR2P)
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

extern "C" void FortranCInterface_GLOBAL(zungqr, ZUNGQR)
  (const int* const M,
   const int* const N,
   const int* const K,
   std::complex<double> A[],
   const int* const LDA,
   std::complex<double> TAU[],
   std::complex<double> WORK[],
   const int* const LWORK,
   int* const INFO);

extern "C" void FortranCInterface_GLOBAL(zunm2r, ZUNM2R)
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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  // If _GEQRFP is available, LAPACK::GEQRF() calls it.  If _LARFP is
  // available, LAPACK::GEQRF() calls _GEQRF, which uses _LARFP.
#ifdef HAVE_LAPACK_ZGEQRFP
  template <>
  const bool LAPACK<int, std::complex<double> >::QR_produces_R_factor_with_nonnegative_diagonal = true;
#else
#  ifdef HAVE_LAPACK_ZLARFP
  template <>
  const bool LAPACK<int, std::complex<double> >::QR_produces_R_factor_with_nonnegative_diagonal = true;
#  else
  template <>
  const bool LAPACK<int, std::complex<double> >::QR_produces_R_factor_with_nonnegative_diagonal = false;
#  endif
#endif

  ////////////////////////////////////////////////////////////////////////////
  // LARFP (implemented with _LARFGP if available, else with _LARFP if
  // available, else fall back to _LARFG)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void 
  LAPACK<int, std::complex<double> >::LARFP (const int n, 
					     std::complex<double>& alpha, 
					     std::complex<double> x[], 
					     const int incx, 
					     std::complex<double>& tau)
  {
#ifdef HAVE_LAPACK_ZLARFGP
    FortranCInterface_GLOBAL(zlarfgp,ZLARFGP) (&n, &alpha, x, &incx, &tau);
#  ifdef HAVE_LAPACK_ZLARFP
    FortranCInterface_GLOBAL(zlarfp,ZLARFP) (&n, &alpha, x, &incx, &tau);
#  else
    FortranCInterface_GLOBAL(zlarfg,ZLARFG) (&n, &alpha, x, &incx, &tau);
#  endif // HAVE_LAPACK_ZLARFP
#endif // HAVE_LAPACK_ZLARFGP
  }

  ////////////////////////////////////////////////////////////////////////////
  // GEQRF (implemented with _GEQRFP if available, else fall back to _GEQRF)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void
  LAPACK<int, std::complex<double> >::GEQRF (const int m,
					     const int n, 
					     std::complex<double> A[],
					     const int lda, 
					     std::complex<double> tau[],
					     std::complex<double> work[],
					     const int lwork,
					     int* const INFO)
  {
#ifdef HAVE_LAPACK_ZGEQRFP
    FortranCInterface_GLOBAL(zgeqrfp, ZGEQRFP) (&m, &n, A, &lda, tau, work, &lwork, INFO);
#else
    FortranCInterface_GLOBAL(zgeqrf, ZGEQRF) (&m, &n, A, &lda, tau, work, &lwork, INFO);
#endif // HAVE_LAPACK_ZGEQRFP
  }

  ////////////////////////////////////////////////////////////////////////////
  // GEQR2 (implemented with _GEQR2P if available, else fall back to _GEQR2)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void
  LAPACK<int, std::complex<double> >::GEQR2 (const int m,
					     const int n, 
					     std::complex<double> A[],
					     const int lda, 
					     std::complex<double> tau[],
					     std::complex<double> work[],
					     int* const INFO)
  {
#ifdef HAVE_LAPACK_ZGEQR2P
    FortranCInterface_GLOBAL(zgeqr2p, ZGEQR2P) (&m, &n, A, &lda, tau, work, INFO);
#else
    FortranCInterface_GLOBAL(zgeqr2, ZGEQR2) (&m, &n, A, &lda, tau, work, INFO);
#endif // HAVE_LAPACK_ZGEQR2P
  }

  template <>
  void
  LAPACK<int, std::complex<double> >::ORM2R (const char* const side,
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
    FortranCInterface_GLOBAL(zunm2r, ZUNM2R) (side, trans, &m, &n, &k, A, &lda, tau, C, &ldc, work, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<double> >::ORGQR (const int m,
					     const int n,
					     const int k,
					     std::complex<double> A[],
					     const int lda,
					     std::complex<double> tau[],
					     std::complex<double> work[],
					     const int lwork,
					     int* const INFO)
  {
    FortranCInterface_GLOBAL(zungqr, ZUNGQR) (&m, &n, &k, A, &lda, tau, work, &lwork, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<double> >::POTRF (const char* const uplo,
					     const int n,
					     std::complex<double> A[],
					     const int lda,
					     int* const INFO)
  {
    FortranCInterface_GLOBAL(zpotrf, ZPOTRF) (uplo, &n, A, &lda, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<double> >::POTRS (const char* const uplo,
					     const int n,
					     const int nrhs,
					     const std::complex<double> A[],
					     const int lda,
					     std::complex<double> B[],
					     const int ldb,
					     int* const INFO)
  {
    FortranCInterface_GLOBAL(zpotrs, ZPOTRS) (uplo, &n, &nrhs, A, &lda, B, &ldb, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<double> >::POTRI (const char* const uplo, 
					     const int n, 
					     std::complex<double> A[], 
					     const int lda, 
					     int* const INFO)
  {
    FortranCInterface_GLOBAL(zpotri, ZPOTRI) (uplo, &n, A, &lda, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<double> >::LARNV (const int idist, 
					     int iseed[],
					     const int n,
					     std::complex<double> x[])
  {
    FortranCInterface_GLOBAL(zlarnv, ZLARNV) (&idist, iseed, &n, x);
  }

} // namespace TSQR
