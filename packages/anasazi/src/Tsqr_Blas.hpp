#ifndef __TSQR_TsqrBlas_hpp
#define __TSQR_TsqrBlas_hpp

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  template< class Ordinal, class Scalar >
  class BLAS {
  public:
    BLAS () {}

    void 
    GEMV (const char* const trans, 
	  const Ordinal m, 
	  const Ordinal n,
	  const Scalar alpha,
	  const Scalar A[],
	  const Ordinal lda,
	  const Scalar x[],
	  const Ordinal incx,
	  const Scalar beta,
	  Scalar y[],
	  const Ordinal incy);

    void
    GEMM (const char* const transa,
	  const char* const transb,
	  const Ordinal m,
	  const Ordinal n,
	  const Ordinal k,
	  const Scalar alpha,
	  const Scalar A[],
	  const Ordinal lda,
	  const Scalar B[],
	  const Ordinal ldb,
	  const Scalar beta,
	  Scalar C[],
	  const Ordinal ldc);

    /// DGER, SGER for real, ZGERC, CGERC for complex (FIXME mfh 28
    /// Apr 2010 is that right?  that's what ZLARF uses...).
    void
    GER (const Ordinal m,
	 const Ordinal n,
	 const Scalar alpha,
	 const Scalar x[],
	 const Ordinal incx,
	 const Scalar y[],
	 const Ordinal incy,
	 Scalar A[],
	 const Ordinal lda);

    void
    TRSM (const char* const side,
	  const char* const uplo,
	  const char* const transa,
	  const char* const diag,
	  const Ordinal m,
	  const Ordinal n,
	  const Scalar alpha,
	  const Scalar A[],
	  const Ordinal lda,
	  Scalar B[],
	  const Ordinal ldb);

  private:
    BLAS (const BLAS&);
    BLAS& operator= (const BLAS&);
  };

} // namespace TSQR

#endif // __TSQR_TsqrBlas_hpp
