#include <Tsqr_Lapack.hpp>
#include <Tsqr_FortranCInterface.hpp>
#include <Tsqr_Config.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

extern "C" void FortranCInterface_GLOBAL(dlarnv, DLARNV)
  (const int* const IDIST,
   int ISEED[],
   const int* const N,
   double X[]);

extern "C" void FortranCInterface_GLOBAL(dpotri, DPOTRI)
  (const char* const UPLO,
   const int* const N,
   double A[],
   const int* const LDA,
   int* const INFO);

extern "C" void FortranCInterface_GLOBAL(dpotrf, DPOTRF)
  (const char* const UPLO,
   const int* const N,
   double A[],
   const int* const LDA,
   int* const INFO);

extern "C" void FortranCInterface_GLOBAL(dpotrs, DPOTRS)
  (const char* const UPLO,
   const int* const N,
   const int* const NRHS,
   const double A[],
   const int* const LDA,
   double B[],
   const int* const LDB,
   int* const INFO);

#ifdef HAVE_LAPACK_DLARFGP
extern "C" void FortranCInterface_GLOBAL(dlarfgp,DLARFGP)
  (const int* const N,    // IN
   double* const ALPHA,   // IN/OUT
   double X[],            // IN/OUT
   const int* const INCX, // IN
   double* const TAU);    // OUT
#else
#  ifdef HAVE_LAPACK_DLARFP
extern "C" void FortranCInterface_GLOBAL(dlarfp,DLARFP)
  (const int* const N,    // IN
   double* const ALPHA,   // IN/OUT
   double X[],            // IN/OUT
   const int* const INCX, // IN
   double* const TAU);    // OUT
#  else
extern "C" void FortranCInterface_GLOBAL(dlarfg,DLARFG)
  (const int* const N,    // IN
   double* const ALPHA,   // IN/OUT
   double X[],            // IN/OUT
   const int* const INCX, // IN
   double* const TAU);    // OUT
#  endif // HAVE_LAPACK_DLARFP
#endif // HAVE_LAPACK_DLARFGP

extern "C" void FortranCInterface_GLOBAL(dgeqrf, DGEQRF)
  (const int* const M,
   const int* const N,
   double A[],
   const int* const LDA,
   double TAU[],
   double WORK[],
   const int* const LWORK,
   int* const INFO);

#ifdef HAVE_LAPACK_DGEQRFP
extern "C" void FortranCInterface_GLOBAL(dgeqrfp, DGEQRFP)
  (const int* const M,
   const int* const N,
   double A[],
   const int* const LDA,
   double TAU[],
   double WORK[],
   const int* const LWORK,
   int* const INFO);
#endif // HAVE_LAPACK_DGEQRFP

extern "C" void FortranCInterface_GLOBAL(dgeqr2, DGEQR2)
  (const int* const M,
   const int* const N,
   double A[],
   const int* const LDA,
   double TAU[],
   double WORK[],
   int* const INFO);

#ifdef HAVE_LAPACK_DGEQR2P
extern "C" void FortranCInterface_GLOBAL(dgeqr2p, DGEQR2P)
  (const int* const M,
   const int* const N,
   double A[],
   const int* const LDA,
   double TAU[],
   double WORK[],
   int* const INFO);
#endif // HAVE_LAPACK_DGEQR2P

extern "C" void FortranCInterface_GLOBAL(dorm2r, DORM2R)
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

extern "C" void FortranCInterface_GLOBAL(dorgqr, DORGQR)
  (const int* const M,
   const int* const N,
   const int* const K,
   double A[],
   const int* const LDA,
   double TAU[],
   double WORK[],
   const int* const LWORK,
   int* const INFO);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  // If _GEQRFP is available, LAPACK::GEQRF() calls it.  If _LARFP is
  // available, LAPACK::GEQRF() calls _GEQRF, which uses _LARFP.
#ifdef HAVE_LAPACK_DGEQRFP
  template <>
  const bool LAPACK<int, double >::QR_produces_R_factor_with_nonnegative_diagonal = true;
#else
#  ifdef HAVE_LAPACK_DLARFP
  template <>
  const bool LAPACK<int, double >::QR_produces_R_factor_with_nonnegative_diagonal = true;
#  else
  template <>
  const bool LAPACK<int, double >::QR_produces_R_factor_with_nonnegative_diagonal = false;
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
    FortranCInterface_GLOBAL(dlarfgp,DLARFGP) (&n, &alpha, x, &incx, &tau);
#  ifdef HAVE_LAPACK_DLARFP
    FortranCInterface_GLOBAL(dlarfp,DLARFP) (&n, &alpha, x, &incx, &tau);
#  else
    FortranCInterface_GLOBAL(dlarfg,DLARFG) (&n, &alpha, x, &incx, &tau);
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
    FortranCInterface_GLOBAL(dgeqrfp, DGEQRFP) (&m, &n, A, &lda, tau, work, &lwork, INFO);
#else
    FortranCInterface_GLOBAL(dgeqrf, DGEQRF) (&m, &n, A, &lda, tau, work, &lwork, INFO);
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
    FortranCInterface_GLOBAL(dgeqr2p, DGEQR2P) (&m, &n, A, &lda, tau, work, INFO);
#else
    FortranCInterface_GLOBAL(dgeqr2, DGEQR2) (&m, &n, A, &lda, tau, work, INFO);
#endif // HAVE_LAPACK_DGEQR2P
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
    FortranCInterface_GLOBAL(dorm2r, DORM2R) (side, trans, &m, &n, &k, A, &lda, tau, C, &ldc, work, INFO);
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
    FortranCInterface_GLOBAL(dorgqr, DORGQR) (&m, &n, &k, A, &lda, tau, work, &lwork, INFO);
  }

  template <>
  void
  LAPACK<int, double >::POTRF (const char* const uplo,
			       const int n,
			       double A[],
			       const int lda,
			       int* const INFO)
  {
    FortranCInterface_GLOBAL(dpotrf, DPOTRF) (uplo, &n, A, &lda, INFO);
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
    FortranCInterface_GLOBAL(dpotrs, DPOTRS) (uplo, &n, &nrhs, A, &lda, B, &ldb, INFO);
  }

  template <>
  void
  LAPACK<int, double >::POTRI (const char* const uplo, 
			       const int n, 
			       double A[], 
			       const int lda, 
			       int* const INFO)
  {
    FortranCInterface_GLOBAL(dpotri, DPOTRI) (uplo, &n, A, &lda, INFO);
  }

  template <>
  void
  LAPACK<int, double >::LARNV (const int idist, 
			       int iseed[],
			       const int n,
			       double x[])
  {
    FortranCInterface_GLOBAL(dlarnv, DLARNV) (&idist, iseed, &n, x);
  }

} // namespace TSQR
