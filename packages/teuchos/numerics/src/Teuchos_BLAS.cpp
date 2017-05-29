// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
// @HEADER

#include "Teuchos_BLAS.hpp"
#include "Teuchos_BLAS_wrappers.hpp"

/* for INTEL_CXML, the second arg may need to be changed to 'one'.  If so
the appropriate declaration of one will need to be added back into
functions that include the macro:
*/

namespace {
#if defined (INTEL_CXML)
        unsigned int one=1;
#endif
} // namespace

#ifdef CHAR_MACRO
#undef CHAR_MACRO
#endif
#if defined (INTEL_CXML)
#define CHAR_MACRO(char_var) &char_var, one
#else
#define CHAR_MACRO(char_var) &char_var
#endif


const char Teuchos::ESideChar[] = {'L' , 'R' };
const char Teuchos::ETranspChar[] = {'N' , 'T' , 'C' };
const char Teuchos::EUploChar[] = {'U' , 'L' };
const char Teuchos::EDiagChar[] = {'U' , 'N' };
const char Teuchos::ETypeChar[] = {'G' , 'L', 'U', 'H', 'B', 'Q', 'Z' };
//const char Teuchos::EFactChar[] = {'F', 'N' };
//const char Teuchos::ENormChar[] = {'O', 'I' };
//const char Teuchos::ECompQChar[] = {'N', 'I', 'V' };
//const char Teuchos::EJobChar[] = {'E', 'V', 'B' };
//const char Teuchos::EJobSChar[] = {'E', 'S' };
//const char Teuchos::EJobVSChar[] = {'V', 'N' };
//const char Teuchos::EHowmnyChar[] = {'A', 'S' };
//const char Teuchos::ECMachChar[] = {'E', 'S', 'B', 'P', 'N', 'R', 'M', 'U', 'L', 'O' };
//const char Teuchos::ESortChar[] = {'N', 'S'};


namespace {


template<typename Scalar>
Scalar generic_dot(const int& n, const Scalar* x, const int& incx,
  const Scalar* y, const int& incy)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  Scalar dot = 0.0;
  if (incx==1 && incy==1) {
    for (int i = 0; i < n; ++i)
      dot += (*x++)*ST::conjugate(*y++);
  }
  else {
    if (incx < 0)
      x = x - incx*(n-1);
    if (incy < 0)
      y = y - incy*(n-1);
    for (int i = 0; i < n; ++i, x+=incx, y+=incy)
      dot += (*x)*ST::conjugate(*y);
  }
  return dot;
}


} // namespace


namespace Teuchos {

//Explicitly instantiating these templates for windows due to an issue with
//resolving them when linking dlls.
#ifdef _MSC_VER
#  ifdef HAVE_TEUCHOS_COMPLEX
     template BLAS<long int, std::complex<float> >;
     template BLAS<long int, std::complex<double> >;
#  endif
     template BLAS<long int, float>;
     template BLAS<long int, double>;
#endif

  // *************************** BLAS<int,float> DEFINITIONS ******************************

  void BLAS<int, float>::ROTG(float* da, float* db, float* c, float* s) const
  { SROTG_F77(da, db, c, s ); }

  void BLAS<int, float>::ROT(const int& n, float* dx, const int& incx, float* dy, const int& incy, float* c, float* s) const
  { SROT_F77(&n, dx, &incx, dy, &incy, c, s); }


  float BLAS<int, float>::ASUM(const int& n, const float* x, const int& incx) const
  {
#if defined(HAVE_TEUCHOS_BLASFLOAT_APPLE_VECLIB_BUGFIX)
    return cblas_sasum(n, x, incx);
#elif defined(HAVE_TEUCHOS_BLASFLOAT)
    float tmp = SASUM_F77(&n, x, &incx);
    return tmp;
#else
    typedef ScalarTraits<float> ST;
    float sum = 0.0;
    if (incx == 1) {
      for (int i = 0; i < n; ++i)
        sum += ST::magnitude(*x++);
    }
    else {
      for (int i = 0; i < n; ++i, x+=incx)
        sum += ST::magnitude(*x);
    }
    return sum;
#endif
  }

  void BLAS<int, float>::AXPY(const int& n, const float& alpha, const float* x, const int& incx, float* y, const int& incy) const
  { SAXPY_F77(&n, &alpha, x, &incx, y, &incy); }

  void BLAS<int, float>::COPY(const int& n, const float* x, const int& incx, float* y, const int& incy) const
  { SCOPY_F77(&n, x, &incx, y, &incy); }

  float BLAS<int, float>::DOT(const int& n, const float* x, const int& incx, const float* y, const int& incy) const
  {
#if defined(HAVE_TEUCHOS_BLASFLOAT_APPLE_VECLIB_BUGFIX)
    return cblas_sdot(n, x, incx, y, incy);
#elif defined(HAVE_TEUCHOS_BLASFLOAT)
    return SDOT_F77(&n, x, &incx, y, &incy);
#else
    return generic_dot(n, x, incx, y, incy);
#endif
  }

  int BLAS<int, float>::IAMAX(const int& n, const float* x, const int& incx) const
  { return ISAMAX_F77(&n, x, &incx); }

  float BLAS<int, float>::NRM2(const int& n, const float* x, const int& incx) const
  {
#if defined(HAVE_TEUCHOS_BLASFLOAT_APPLE_VECLIB_BUGFIX)
    return cblas_snrm2(n, x, incx);
#elif defined(HAVE_TEUCHOS_BLASFLOAT)
    return SNRM2_F77(&n, x, &incx);
#else
    return ScalarTraits<float>::squareroot(generic_dot(n, x, incx, x, incx));
#endif
  }

  void BLAS<int, float>::SCAL(const int& n, const float& alpha, float* x, const int& incx) const
  { SSCAL_F77(&n, &alpha, x, &incx); }

  void BLAS<int, float>::GEMV(ETransp trans, const int& m, const int& n, const float& alpha, const float* A, const int& lda, const float* x, const int& incx, const float& beta, float* y, const int& incy) const
  { SGEMV_F77(CHAR_MACRO(ETranspChar[trans]), &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy); }

  void BLAS<int, float>::GER(const int& m, const int& n, const float& alpha, const float* x, const int& incx, const float* y, const int& incy, float* A, const int& lda) const
  { SGER_F77(&m, &n, &alpha, x, &incx, y, &incy, A, &lda); }

  void BLAS<int, float>::TRMV(EUplo uplo, ETransp trans, EDiag diag, const int& n, const float* A, const int& lda, float* x, const int& incx) const
  { STRMV_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), CHAR_MACRO(EDiagChar[diag]), &n, A, &lda, x, &incx); }

  void BLAS<int, float>::GEMM(ETransp transa, ETransp transb, const int& m, const int& n, const int& k, const float& alpha, const float* A, const int& lda, const float* B, const int& ldb, const float& beta, float* C, const int& ldc) const
  { SGEMM_F77(CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(ETranspChar[transb]), &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }

  void BLAS<int, float>::SWAP(const int& n, float* const x, const int& incx, float* const y, const int& incy) const
  {
    SSWAP_F77 (&n, x, &incx, y, &incy);
  }

  void BLAS<int, float>::SYMM(ESide side, EUplo uplo, const int& m, const int& n, const float& alpha, const float* A, const int& lda, const float* B, const int& ldb, const float& beta, float* C, const int& ldc) const
  { SSYMM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }

  void BLAS<int, float>::SYRK(EUplo uplo, ETransp trans, const int& n, const int& k, const float& alpha, const float* A, const int& lda, const float& beta, float* C, const int& ldc) const
  { SSYRK_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), &n, &k, &alpha, A, &lda, &beta, C, &ldc); }

  void BLAS<int, float>::HERK(EUplo uplo, ETransp trans, const int& n, const int& k, const float& alpha, const float* A, const int& lda, const float& beta, float* C, const int& ldc) const
  { SSYRK_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), &n, &k, &alpha, A, &lda, &beta, C, &ldc); }

  void BLAS<int, float>::TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int& m, const int& n, const float& alpha, const float* A, const int& lda, float* B, const int& ldb) const
  { STRMM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(EDiagChar[diag]), &m, &n, &alpha, A, &lda, B, &ldb); }

  void BLAS<int, float>::TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int& m, const int& n, const float& alpha, const float* A, const int& lda, float* B, const int& ldb) const
  { STRSM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(EDiagChar[diag]), &m, &n, &alpha, A, &lda, B, &ldb); }

  // *************************** BLAS<int,double> DEFINITIONS ******************************

  void BLAS<int, double>::ROTG(double* da, double* db, double* c, double* s) const
  { DROTG_F77(da, db, c, s); }

  void BLAS<int, double>::ROT(const int& n, double* dx, const int& incx, double* dy, const int& incy, double* c, double* s) const
  { DROT_F77(&n, dx, &incx, dy, &incy, c, s); }

  double BLAS<int, double>::ASUM(const int& n, const double* x, const int& incx) const
  { return DASUM_F77(&n, x, &incx); }

  void BLAS<int, double>::AXPY(const int& n, const double& alpha, const double* x, const int& incx, double* y, const int& incy) const
  { DAXPY_F77(&n, &alpha, x, &incx, y, &incy); }

  void BLAS<int, double>::COPY(const int& n, const double* x, const int& incx, double* y, const int& incy) const
  { DCOPY_F77(&n, x, &incx, y, &incy); }

  double BLAS<int, double>::DOT(const int& n, const double* x, const int& incx, const double* y, const int& incy) const
  {
    return DDOT_F77(&n, x, &incx, y, &incy);
  }

  int BLAS<int, double>::IAMAX(const int& n, const double* x, const int& incx) const
  { return IDAMAX_F77(&n, x, &incx); }

  double BLAS<int, double>::NRM2(const int& n, const double* x, const int& incx) const
  { return DNRM2_F77(&n, x, &incx); }

  void BLAS<int, double>::SCAL(const int& n, const double& alpha, double* x, const int& incx) const
  { DSCAL_F77(&n, &alpha, x, &incx); }

  void BLAS<int, double>::GEMV(ETransp trans, const int& m, const int& n, const double& alpha, const double* A, const int& lda, const double* x, const int& incx, const double& beta, double* y, const int& incy) const
  { DGEMV_F77(CHAR_MACRO(ETranspChar[trans]), &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy); }

  void BLAS<int, double>::GER(const int& m, const int& n, const double& alpha, const double* x, const int& incx, const double* y, const int& incy, double* A, const int& lda) const
  { DGER_F77(&m, &n, &alpha, x, &incx, y, &incy, A, &lda); }

  void BLAS<int, double>::TRMV(EUplo uplo, ETransp trans, EDiag diag, const int& n, const double* A, const int& lda, double* x, const int& incx) const
  { DTRMV_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), CHAR_MACRO(EDiagChar[diag]), &n, A, &lda, x, &incx); }

  void BLAS<int, double>::GEMM(ETransp transa, ETransp transb, const int& m, const int& n, const int& k, const double& alpha, const double* A, const int& lda, const double* B, const int& ldb, const double& beta, double* C, const int& ldc) const
  { DGEMM_F77(CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(ETranspChar[transb]), &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }

  void BLAS<int, double>::SWAP(const int& n, double* const x, const int& incx, double* const y, const int& incy) const
  {
    DSWAP_F77 (&n, x, &incx, y, &incy);
  }

  void BLAS<int, double>::SYMM(ESide side, EUplo uplo, const int& m, const int& n, const double& alpha, const double* A, const int& lda, const double* B, const int& ldb, const double& beta, double* C, const int& ldc) const
  { DSYMM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }

  void BLAS<int, double>::SYRK(EUplo uplo, ETransp trans, const int& n, const int& k, const double& alpha, const double* A, const int& lda, const double& beta, double* C, const int& ldc) const
  { DSYRK_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), &n, &k, &alpha, A, &lda, &beta, C, &ldc); }

  void BLAS<int, double>::HERK(EUplo uplo, ETransp trans, const int& n, const int& k, const double& alpha, const double* A, const int& lda, const double& beta, double* C, const int& ldc) const
  { DSYRK_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), &n, &k, &alpha, A, &lda, &beta, C, &ldc); }

  void BLAS<int, double>::TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int& m, const int& n, const double& alpha, const double* A, const int& lda, double* B, const int& ldb) const
  { DTRMM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(EDiagChar[diag]), &m, &n, &alpha, A, &lda, B, &ldb); }

  void BLAS<int, double>::TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int& m, const int& n, const double& alpha, const double* A, const int& lda, double* B, const int& ldb) const
  { DTRSM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(EDiagChar[diag]), &m, &n, &alpha, A, &lda, B, &ldb); }

#ifdef HAVE_TEUCHOS_COMPLEX

  // *************************** BLAS<int,std::complex<float> > DEFINITIONS ******************************

  void BLAS<int, std::complex<float> >::ROTG(std::complex<float>* da, std::complex<float>* db, float* c, std::complex<float>* s) const
  { CROTG_F77(da, db, c, s ); }

  void BLAS<int, std::complex<float> >::ROT(const int& n, std::complex<float>* dx, const int& incx, std::complex<float>* dy, const int& incy, float* c, std::complex<float>* s) const
  { CROT_F77(&n, dx, &incx, dy, &incy, c, s); }

  float BLAS<int, std::complex<float> >::ASUM(const int& n, const std::complex<float>* x, const int& incx) const
  {
#if defined(HAVE_TEUCHOS_BLASFLOAT_APPLE_VECLIB_BUGFIX)
    return cblas_scasum(n, x, incx);
#elif defined(HAVE_TEUCHOS_BLASFLOAT_DOUBLE_RETURN)
    return (float) SCASUM_F77(&n, x, &incx);
#elif defined(HAVE_TEUCHOS_BLASFLOAT)
    return SCASUM_F77(&n, x, &incx);
#else // Wow, you just plain don't have this routine.
    // mfh 01 Feb 2013: See www.netlib.org/blas/scasum.f.
    // I've enhanced this by accumulating in double precision.
    double result = 0;
    if (incx == 1) {
      for (int i = 0; i < n; ++i) {
        result += std::abs (std::real (x[i])) + std::abs (std::imag (x[i]));
      }
    } else {
      const int nincx = n * incx;
      for (int i = 0; i < nincx; i += incx) {
        result += std::abs (std::real (x[i])) + std::abs (std::imag (x[i]));
      }
    }
    return static_cast<float> (result);
#endif
  }

  void BLAS<int, std::complex<float> >::AXPY(const int& n, const std::complex<float> alpha, const std::complex<float>* x, const int& incx, std::complex<float>* y, const int& incy) const
  { CAXPY_F77(&n, &alpha, x, &incx, y, &incy); }

  void BLAS<int, std::complex<float> >::COPY(const int& n, const std::complex<float>* x, const int& incx, std::complex<float>* y, const int& incy) const
  { CCOPY_F77(&n, x, &incx, y, &incy); }

  std::complex<float> BLAS<int, std::complex<float> >::DOT(const int& n, const std::complex<float>* x, const int& incx, const std::complex<float>* y, const int& incy) const
  {
#if defined(HAVE_TEUCHOS_BLASFLOAT_APPLE_VECLIB_BUGFIX)
    std::complex<float> z;
    cblas_cdotc_sub(n,x,incx,y,incy,&z);
    return z;
#elif defined(HAVE_COMPLEX_BLAS_PROBLEM) && defined(HAVE_FIXABLE_COMPLEX_BLAS_PROBLEM)
    std::complex<float> z;
    CDOT_F77(&z, &n, x, &incx, y, &incy);
    return z;
#elif defined(HAVE_TEUCHOS_BLASFLOAT)
    return CDOT_F77(&n, x, &incx, y, &incy);
#else // Wow, you just plain don't have this routine.
    // mfh 01 Feb 2013: See www.netlib.org/blas/cdotc.f.
    // I've enhanced this by accumulating in double precision.
    std::complex<double> result (0, 0);
    if (n >= 0) {
      if (incx == 1 && incy == 1) {
        for (int i = 0; i < n; ++i) {
          result += std::conj (x[i]) * y[i];
        }
      } else {
        int ix = 0;
        int iy = 0;
        if (incx < 0) {
          ix = (1-n) * incx;
        }
        if (incy < 0) {
          iy = (1-n) * incy;
        }
        for (int i = 0; i < n; ++i) {
          result += std::conj (x[ix]) * y[iy];
          ix += incx;
          iy += incy;
        }
      }
    }
    return static_cast<std::complex<float> > (result);
#endif
  }

  int BLAS<int, std::complex<float> >::IAMAX(const int& n, const std::complex<float>* x, const int& incx) const
  { return ICAMAX_F77(&n, x, &incx); }

  float BLAS<int, std::complex<float> >::NRM2(const int& n, const std::complex<float>* x, const int& incx) const
  {
#if defined(HAVE_TEUCHOS_BLASFLOAT_APPLE_VECLIB_BUGFIX)
    return cblas_scnrm2(n, x, incx);
#elif defined(HAVE_TEUCHOS_BLASFLOAT_DOUBLE_RETURN)
    return (float) SCNRM2_F77(&n, x, &incx);
#elif defined(HAVE_TEUCHOS_BLASFLOAT)
    return SCNRM2_F77(&n, x, &incx);
#else // Wow, you just plain don't have this routine.
    // mfh 01 Feb 2013: See www.netlib.org/blas/scnrm2.f.
    // I've enhanced this by accumulating in double precision.
    if (n < 1 || incx < 1) {
      return 0;
    } else {
      double scale = 0;
      double ssq = 1;

      const int upper = 1 + (n-1)*incx;
      for (int ix = 0; ix < upper; ix += incx) {
        // The reference BLAS implementation cleverly scales the
        // intermediate result. so that even if the square of the norm
        // would overflow, computing the norm itself does not.  Hence,
        // "ssq" for "scaled square root."
        if (std::real (x[ix]) != 0) {
          const double temp = std::abs (std::real (x[ix]));
          if (scale < temp) {
            const double scale_over_temp = scale / temp;
            ssq = 1 + ssq * scale_over_temp*scale_over_temp;
            // New scaling factor: biggest (in magnitude) real or imaginary part seen thus far.
            scale = temp;
          } else {
            const double temp_over_scale = temp / scale;
            ssq = ssq + temp_over_scale*temp_over_scale;
          }
        }
        if (std::imag (x[ix]) != 0) {
          const double temp = std::abs (std::imag (x[ix]));
          if (scale < temp) {
            const double scale_over_temp = scale / temp;
            ssq = 1 + ssq * scale_over_temp*scale_over_temp;
            // New scaling factor: biggest (in magnitude) real or imaginary part seen thus far.
            scale = temp;
          } else {
            const double temp_over_scale = temp / scale;
            ssq = ssq + temp_over_scale*temp_over_scale;
          }
        }
      }
      return static_cast<float> (scale * std::sqrt (ssq));
    }
#endif
  }

  void BLAS<int, std::complex<float> >::SCAL(const int& n, const std::complex<float> alpha, std::complex<float>* x, const int& incx) const
  { CSCAL_F77(&n, &alpha, x, &incx); }

  void BLAS<int, std::complex<float> >::GEMV(ETransp trans, const int& m, const int& n, const std::complex<float> alpha, const std::complex<float>* A, const int& lda, const std::complex<float>* x, const int& incx, const std::complex<float> beta, std::complex<float>* y, const int& incy) const
  { CGEMV_F77(CHAR_MACRO(ETranspChar[trans]), &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy); }

  void BLAS<int, std::complex<float> >::GER(const int& m, const int& n, const std::complex<float> alpha, const std::complex<float>* x, const int& incx, const std::complex<float>* y, const int& incy, std::complex<float>* A, const int& lda) const
  { CGER_F77(&m, &n, &alpha, x, &incx, y, &incy, A, &lda); }

  void BLAS<int, std::complex<float> >::TRMV(EUplo uplo, ETransp trans, EDiag diag, const int& n, const std::complex<float>* A, const int& lda, std::complex<float>* x, const int& incx) const
  { CTRMV_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), CHAR_MACRO(EDiagChar[diag]), &n, A, &lda, x, &incx); }

  void BLAS<int, std::complex<float> >::GEMM(ETransp transa, ETransp transb, const int& m, const int& n, const int& k, const std::complex<float> alpha, const std::complex<float>* A, const int& lda, const std::complex<float>* B, const int& ldb, const std::complex<float> beta, std::complex<float>* C, const int& ldc) const
  { CGEMM_F77(CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(ETranspChar[transb]), &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }

  void BLAS<int, std::complex<float> >::SWAP(const int& n, std::complex<float>* const x, const int& incx, std::complex<float>* const y, const int& incy) const
  {
    CSWAP_F77 (&n, x, &incx, y, &incy);
  }

  void BLAS<int, std::complex<float> >::SYMM(ESide side, EUplo uplo, const int& m, const int& n, const std::complex<float> alpha, const std::complex<float>* A, const int& lda, const std::complex<float>* B, const int& ldb, const std::complex<float> beta, std::complex<float>* C, const int& ldc) const
  { CSYMM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }

  void BLAS<int, std::complex<float> >::SYRK(EUplo uplo, ETransp trans, const int& n, const int& k, const std::complex<float> alpha, const std::complex<float>* A, const int& lda, const std::complex<float> beta, std::complex<float>* C, const int& ldc) const
  { CSYRK_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), &n, &k, &alpha, A, &lda, &beta, C, &ldc); }

  void BLAS<int, std::complex<float> >::HERK(EUplo uplo, ETransp trans, const int& n, const int& k, const std::complex<float> alpha, const std::complex<float>* A, const int& lda, const std::complex<float> beta, std::complex<float>* C, const int& ldc) const
  { CHERK_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), &n, &k, &alpha, A, &lda, &beta, C, &ldc); }

  void BLAS<int, std::complex<float> >::TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int& m, const int& n, const std::complex<float> alpha, const std::complex<float>* A, const int& lda, std::complex<float>* B, const int& ldb) const
  { CTRMM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(EDiagChar[diag]), &m, &n, &alpha, A, &lda, B, &ldb); }

  void BLAS<int, std::complex<float> >::TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int& m, const int& n, const std::complex<float> alpha, const std::complex<float>* A, const int& lda, std::complex<float>* B, const int& ldb) const
  { CTRSM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(EDiagChar[diag]), &m, &n, &alpha, A, &lda, B, &ldb); }

  // *************************** BLAS<int,std::complex<double> > DEFINITIONS ******************************

  void BLAS<int, std::complex<double> >::ROTG(std::complex<double>* da, std::complex<double>* db, double* c, std::complex<double>* s) const
  { ZROTG_F77(da, db, c, s); }

  void BLAS<int, std::complex<double> >::ROT(const int& n, std::complex<double>* dx, const int& incx, std::complex<double>* dy, const int& incy, double* c, std::complex<double>* s) const
  { ZROT_F77(&n, dx, &incx, dy, &incy, c, s); }

  double BLAS<int, std::complex<double> >::ASUM(const int& n, const std::complex<double>* x, const int& incx) const
  { return ZASUM_F77(&n, x, &incx); }

  void BLAS<int, std::complex<double> >::AXPY(const int& n, const std::complex<double> alpha, const std::complex<double>* x, const int& incx, std::complex<double>* y, const int& incy) const
  { ZAXPY_F77(&n, &alpha, x, &incx, y, &incy); }

  void BLAS<int, std::complex<double> >::COPY(const int& n, const std::complex<double>* x, const int& incx, std::complex<double>* y, const int& incy) const
  { ZCOPY_F77(&n, x, &incx, y, &incy); }

  std::complex<double> BLAS<int, std::complex<double> >::DOT(const int& n, const std::complex<double>* x, const int& incx, const std::complex<double>* y, const int& incy) const
  {
#if defined(HAVE_TEUCHOS_BLASFLOAT_APPLE_VECLIB_BUGFIX)
    std::complex<double> z;
    cblas_zdotc_sub(n,x,incx,y,incy,&z);
    return z;
#elif defined(HAVE_COMPLEX_BLAS_PROBLEM)
#  if defined(HAVE_FIXABLE_COMPLEX_BLAS_PROBLEM)
    std::complex<double> z;
    ZDOT_F77(&z, &n, x, &incx, y, &incy);
    return z;
#  else
    // mfh 01 Feb 2013: Your complex BLAS is broken, but the problem
    // doesn't have the easy workaround.  I'll just reimplement the
    // missing routine here.  See www.netlib.org/blas/zdotc.f.
    std::complex<double> ztemp (0, 0);
    if (n > 0) {
      if (incx == 1 && incy == 1) {
        for (int i = 0; i < n; ++i) {
          ztemp += std::conj (x[i]) * y[i];
        }
      } else {
        int ix = 0;
        int iy = 0;
        if (incx < 0) {
          ix = (1-n)*incx;
        }
        if (incy < 0) {
          iy = (1-n)*incy;
        }
        for (int i = 0; i < n; ++i) {
          ztemp += std::conj (x[ix]) * y[iy];
          ix += incx;
          iy += incy;
        }
      }
    }
    return ztemp;

#  endif // defined(HAVE_FIXABLE_COMPLEX_BLAS_PROBLEM)
#else
    return ZDOT_F77(&n, x, &incx, y, &incy);
#endif
  }

  int BLAS<int, std::complex<double> >::IAMAX(const int& n, const std::complex<double>* x, const int& incx) const
  { return IZAMAX_F77(&n, x, &incx); }

  double BLAS<int, std::complex<double> >::NRM2(const int& n, const std::complex<double>* x, const int& incx) const
  { return ZNRM2_F77(&n, x, &incx); }

  void BLAS<int, std::complex<double> >::SCAL(const int& n, const std::complex<double> alpha, std::complex<double>* x, const int& incx) const
  { ZSCAL_F77(&n, &alpha, x, &incx); }

  void BLAS<int, std::complex<double> >::GEMV(ETransp trans, const int& m, const int& n, const std::complex<double> alpha, const std::complex<double>* A, const int& lda, const std::complex<double>* x, const int& incx, const std::complex<double> beta, std::complex<double>* y, const int& incy) const
  { ZGEMV_F77(CHAR_MACRO(ETranspChar[trans]), &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy); }

  void BLAS<int, std::complex<double> >::GER(const int& m, const int& n, const std::complex<double> alpha, const std::complex<double>* x, const int& incx, const std::complex<double>* y, const int& incy, std::complex<double>* A, const int& lda) const
  { ZGER_F77(&m, &n, &alpha, x, &incx, y, &incy, A, &lda); }

  void BLAS<int, std::complex<double> >::TRMV(EUplo uplo, ETransp trans, EDiag diag, const int& n, const std::complex<double>* A, const int& lda, std::complex<double>* x, const int& incx) const
  { ZTRMV_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), CHAR_MACRO(EDiagChar[diag]), &n, A, &lda, x, &incx); }

  void BLAS<int, std::complex<double> >::GEMM(ETransp transa, ETransp transb, const int& m, const int& n, const int& k, const std::complex<double> alpha, const std::complex<double>* A, const int& lda, const std::complex<double>* B, const int& ldb, const std::complex<double> beta, std::complex<double>* C, const int& ldc) const
  { ZGEMM_F77(CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(ETranspChar[transb]), &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }

  void BLAS<int, std::complex<double> >::SWAP(const int& n, std::complex<double>* const x, const int& incx, std::complex<double>* const y, const int& incy) const
  {
    ZSWAP_F77 (&n, x, &incx, y, &incy);
  }

  void BLAS<int, std::complex<double> >::SYMM(ESide side, EUplo uplo, const int& m, const int& n, const std::complex<double> alpha, const std::complex<double>* A, const int& lda, const std::complex<double> *B, const int& ldb, const std::complex<double> beta, std::complex<double> *C, const int& ldc) const
  { ZSYMM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }

  void BLAS<int, std::complex<double> >::SYRK(EUplo uplo, ETransp trans, const int& n, const int& k, const std::complex<double> alpha, const std::complex<double>* A, const int& lda, const std::complex<double> beta, std::complex<double>* C, const int& ldc) const
  { ZSYRK_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), &n, &k, &alpha, A, &lda, &beta, C, &ldc); }

  void BLAS<int, std::complex<double> >::HERK(EUplo uplo, ETransp trans, const int& n, const int& k, const std::complex<double> alpha, const std::complex<double>* A, const int& lda, const std::complex<double> beta, std::complex<double>* C, const int& ldc) const
  { ZHERK_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), &n, &k, &alpha, A, &lda, &beta, C, &ldc); }

  void BLAS<int, std::complex<double> >::TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int& m, const int& n, const std::complex<double> alpha, const std::complex<double>* A, const int& lda, std::complex<double>* B, const int& ldb) const
  { ZTRMM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(EDiagChar[diag]), &m, &n, &alpha, A, &lda, B, &ldb); }

  void BLAS<int, std::complex<double> >::TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int& m, const int& n, const std::complex<double> alpha, const std::complex<double>* A, const int& lda, std::complex<double>* B, const int& ldb) const
  { ZTRSM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(EDiagChar[diag]), &m, &n, &alpha, A, &lda, B, &ldb); }

#endif // HAVE_TEUCHOS_COMPLEX

}
