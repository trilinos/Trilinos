#include <Tsqr_Lapack.hpp>
#include <Tsqr_FortranCInterface.hpp>
#include <Tsqr_Config.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

extern "C" void FortranCInterface_GLOBAL(slarnv, SLARNV)
  (const int* const IDIST,
   int ISEED[],
   const int* const N,
   float X[]);

extern "C" void FortranCInterface_GLOBAL(spotri, SPOTRI)
  (const char* const UPLO,
   const int* const N,
   float A[],
   const int* const LDA,
   int* const INFO);

extern "C" void FortranCInterface_GLOBAL(spotrf, SPOTRF)
  (const char* const UPLO,
   const int* const N,
   float A[],
   const int* const LDA,
   int* const INFO);

extern "C" void FortranCInterface_GLOBAL(spotrs, SPOTRS)
  (const char* const UPLO,
   const int* const N,
   const int* const NRHS,
   const float A[],
   const int* const LDA,
   float B[],
   const int* const LDB,
   int* const INFO);

#ifdef HAVE_LAPACK_SLARFGP
extern "C" void FortranCInterface_GLOBAL(slarfgp,SLARFGP)
  (const int* const N,    // IN
   float* const ALPHA,   // IN/OUT
   float X[],            // IN/OUT
   const int* const INCX, // IN
   float* const TAU);    // OUT
#else
#  ifdef HAVE_LAPACK_SLARFP
extern "C" void FortranCInterface_GLOBAL(slarfp,SLARFP)
  (const int* const N,    // IN
   float* const ALPHA,   // IN/OUT
   float X[],            // IN/OUT
   const int* const INCX, // IN
   float* const TAU);    // OUT
#  else
extern "C" void FortranCInterface_GLOBAL(slarfg,SLARFG)
  (const int* const N,    // IN
   float* const ALPHA,   // IN/OUT
   float X[],            // IN/OUT
   const int* const INCX, // IN
   float* const TAU);    // OUT
#  endif // HAVE_LAPACK_SLARFP
#endif // HAVE_LAPACK_SLARFGP

extern "C" void FortranCInterface_GLOBAL(sgeqrf, SGEQRF)
  (const int* const M,
   const int* const N,
   float A[],
   const int* const LDA,
   float TAU[],
   float WORK[],
   const int* const LWORK,
   int* const INFO);

#ifdef HAVE_LAPACK_SGEQRFP
extern "C" void FortranCInterface_GLOBAL(sgeqrfp, SGEQRFP)
  (const int* const M,
   const int* const N,
   float A[],
   const int* const LDA,
   float TAU[],
   float WORK[],
   const int* const LWORK,
   int* const INFO);
#endif // HAVE_LAPACK_SGEQRFP

extern "C" void FortranCInterface_GLOBAL(sgeqr2, SGEQR2)
  (const int* const M,
   const int* const N,
   float A[],
   const int* const LDA,
   float TAU[],
   float WORK[],
   int* const INFO);

#ifdef HAVE_LAPACK_SGEQR2P
extern "C" void FortranCInterface_GLOBAL(sgeqr2p, SGEQR2P)
  (const int* const M,
   const int* const N,
   float A[],
   const int* const LDA,
   float TAU[],
   float WORK[],
   int* const INFO);
#endif // HAVE_LAPACK_SGEQR2P

extern "C" void FortranCInterface_GLOBAL(sorm2r, SORM2R)
  (const char* const SIDE,
   const char* const TRANS,
   const int* const M,
   const int* const N,
   const int* const K,
   const float A[],
   const int* const LDA,
   const float TAU[],
   float C[],
   const int* const LDC,
   float WORK[],
   int* const INFO);

extern "C" void FortranCInterface_GLOBAL(sorgqr, SORGQR)
  (const int* const M,
   const int* const N,
   const int* const K,
   float A[],
   const int* const LDA,
   float TAU[],
   float WORK[],
   const int* const LWORK,
   int* const INFO);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  // If _GEQRFP is available, LAPACK::GEQRF() calls it.  If _LARFP is
  // available, LAPACK::GEQRF() calls _GEQRF, which uses _LARFP.
#ifdef HAVE_LAPACK_SGEQRFP
  template <>
  const bool LAPACK<int, float >::QR_produces_R_factor_with_nonnegative_diagonal = true;
#else
#  ifdef HAVE_LAPACK_SLARFP
  template <>
  const bool LAPACK<int, float >::QR_produces_R_factor_with_nonnegative_diagonal = true;
#  else
  template <>
  const bool LAPACK<int, float >::QR_produces_R_factor_with_nonnegative_diagonal = false;
#  endif
#endif

  ////////////////////////////////////////////////////////////////////////////
  // LARFP (implemented with _LARFGP if available, else with _LARFP if
  // available, else fall back to _LARFG)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void 
  LAPACK<int, float >::LARFP (const int n, 
			      float& alpha, 
			      float x[], 
			      const int incx, 
			      float& tau)
  {
#ifdef HAVE_LAPACK_SLARFGP
    FortranCInterface_GLOBAL(slarfgp,SLARFGP) (&n, &alpha, x, &incx, &tau);
#  ifdef HAVE_LAPACK_SLARFP
    FortranCInterface_GLOBAL(slarfp,SLARFP) (&n, &alpha, x, &incx, &tau);
#  else
    FortranCInterface_GLOBAL(slarfg,SLARFG) (&n, &alpha, x, &incx, &tau);
#  endif // HAVE_LAPACK_SLARFP
#endif // HAVE_LAPACK_SLARFGP
  }

  ////////////////////////////////////////////////////////////////////////////
  // GEQRF (implemented with _GEQRFP if available, else fall back to _GEQRF)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void
  LAPACK<int, float >::GEQRF (const int m,
			      const int n, 
			      float A[],
			      const int lda, 
			      float tau[],
			      float work[],
			      const int lwork,
			      int* const INFO)
  {
#ifdef HAVE_LAPACK_SGEQRFP
    FortranCInterface_GLOBAL(sgeqrfp, SGEQRFP) (&m, &n, A, &lda, tau, work, &lwork, INFO);
#else
    FortranCInterface_GLOBAL(sgeqrf, SGEQRF) (&m, &n, A, &lda, tau, work, &lwork, INFO);
#endif // HAVE_LAPACK_SGEQRFP
  }

  ////////////////////////////////////////////////////////////////////////////
  // GEQR2 (implemented with _GEQR2P if available, else fall back to _GEQR2)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void
  LAPACK<int, float >::GEQR2 (const int m,
			      const int n, 
			      float A[],
			      const int lda, 
			      float tau[],
			      float work[],
			      int* const INFO)
  {
#ifdef HAVE_LAPACK_SGEQR2P
    FortranCInterface_GLOBAL(sgeqr2p, SGEQR2P) (&m, &n, A, &lda, tau, work, INFO);
#else
    FortranCInterface_GLOBAL(sgeqr2, SGEQR2) (&m, &n, A, &lda, tau, work, INFO);
#endif // HAVE_LAPACK_SGEQR2P
  }

  template <>
  void
  LAPACK<int, float >::ORM2R (const char* const side,
			      const char* const trans,
			      const int m,
			      const int n,
			      const int k,
			      const float A[],
			      const int lda,
			      const float tau[],
			      float C[],
			      const int ldc,
			      float work[],
			      int* const INFO)
  {
    FortranCInterface_GLOBAL(sorm2r, SORM2R) (side, trans, &m, &n, &k, A, &lda, tau, C, &ldc, work, INFO);
  }

  template <>
  void
  LAPACK<int, float >::ORGQR (const int m,
			      const int n,
			      const int k,
			      float A[],
			      const int lda,
			      float tau[],
			      float work[],
			      const int lwork,
			      int* const INFO)
  {
    FortranCInterface_GLOBAL(sorgqr, SORGQR) (&m, &n, &k, A, &lda, tau, work, &lwork, INFO);
  }

  template <>
  void
  LAPACK<int, float >::POTRF (const char* const uplo,
			      const int n,
			      float A[],
			      const int lda,
			      int* const INFO)
  {
    FortranCInterface_GLOBAL(spotrf, SPOTRF) (uplo, &n, A, &lda, INFO);
  }

  template <>
  void
  LAPACK<int, float >::POTRS (const char* const uplo,
			      const int n,
			      const int nrhs,
			      const float A[],
			      const int lda,
			      float B[],
			      const int ldb,
			      int* const INFO)
  {
    FortranCInterface_GLOBAL(spotrs, SPOTRS) (uplo, &n, &nrhs, A, &lda, B, &ldb, INFO);
  }

  template <>
  void
  LAPACK<int, float >::POTRI (const char* const uplo, 
			      const int n, 
			      float A[], 
			      const int lda, 
			      int* const INFO)
  {
    FortranCInterface_GLOBAL(spotri, SPOTRI) (uplo, &n, A, &lda, INFO);
  }

  template <>
  void
  LAPACK<int, float >::LARNV (const int idist, 
			      int iseed[],
			      const int n,
			      float x[])
  {
    FortranCInterface_GLOBAL(slarnv, SLARNV) (&idist, iseed, &n, x);
  }

} // namespace TSQR
