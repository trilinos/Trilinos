#ifndef __TSQR_Tsqr_Lapack_hpp
#define __TSQR_Tsqr_Lapack_hpp

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  template< class Ordinal, class Scalar >
  class LAPACK {
  public:
    LAPACK () {}

    /// Whether or not the QR factorizations computed by LAPACK::GEQRF()
    /// and LAPACK::GEQR2() produce an R factor with all nonnegative
    /// diagonal entries.  This also corresponds to whether
    /// LAPACK::LARFP() always produces a nonnegative BETA output, and
    /// therefore whether the QR factorizations in the TSQR::Combine
    /// class produce R factors with all negative diagonal entries.
    static const bool QR_produces_R_factor_with_nonnegative_diagonal;

    void 
    LARFP (const Ordinal n, 
	   Scalar& alpha, 
	   Scalar x[], 
	   const Ordinal incx, 
	   Scalar& tau);

    void
    GEQRF  (const Ordinal m,
	    const Ordinal n, 
	    Scalar A[],
	    const Ordinal lda,
	    Scalar tau[],
	    Scalar work[],
	    const int lwork,
	    int* const INFO);

    void 
    GEQR2 (const Ordinal m, 
	   const Ordinal n, 
	   Scalar A[],
	   const Ordinal lda, 
	   Scalar tau[],
	   Scalar work[],
	   int* const INFO);

    void
    ORM2R (const char* const side,
	   const char* const trans,
	   const Ordinal m,
	   const Ordinal n,
	   const Ordinal k,
	   const Scalar A[],
	   const Ordinal lda,
	   const Scalar tau[],
	   Scalar C[],
	   const Ordinal ldc,
	   Scalar work[],
	   int* const info);

    void
    ORGQR (const Ordinal m,
	   const Ordinal n,
	   const Ordinal k,
	   Scalar A[],
	   const Ordinal lda,
	   Scalar tau[],
	   Scalar work[],
	   const int lwork,
	   int* const INFO);

    void
    POTRF (const char* const uplo,
	   const Ordinal n,
	   Scalar A[],
	   const Ordinal lda,
	   int* const INFO);

    void
    POTRS (const char* const uplo,
	   const Ordinal n,
	   const Ordinal nrhs,
	   const Scalar A[],
	   const Ordinal lda,
	   Scalar B[],
	   const Ordinal ldb,
	   int* const INFO);

    void
    POTRI (const char* const uplo, 
	   const Ordinal n, 
	   Scalar A[], 
	   const Ordinal lda, 
	   int* const INFO);

  private:
    LAPACK (const LAPACK&);
    LAPACK& operator= (const LAPACK&);
  };

} // namespace TSQR

#endif // __TSQR_Tsqr_Lapack_hpp
