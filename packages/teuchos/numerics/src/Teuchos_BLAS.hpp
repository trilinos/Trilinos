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

// Kris
// 06.16.03 -- Start over from scratch
// 06.16.03 -- Initial templatization (Tpetra_BLAS.cpp is no longer needed)
// 06.18.03 -- Changed xxxxx_() function calls to XXXXX_F77()
//          -- Added warning messages for generic calls
// 07.08.03 -- Move into Teuchos package/namespace
// 07.24.03 -- The first iteration of BLAS generics is nearing completion. Caveats:
//             * TRSM isn't finished yet; it works for one case at the moment (left side, upper tri., no transpose, no unit diag.)
//             * Many of the generic implementations are quite inefficient, ugly, or both. I wrote these to be easy to debug, not for efficiency or legibility. The next iteration will improve both of these aspects as much as possible.
//             * Very little verification of input parameters is done, save for the character-type arguments (TRANS, etc.) which is quite robust.
//             * All of the routines that make use of both an incx and incy parameter (which includes much of the L1 BLAS) are set up to work iff incx == incy && incx > 0. Allowing for differing/negative values of incx/incy should be relatively trivial.
//             * All of the L2/L3 routines assume that the entire matrix is being used (that is, if A is mxn, lda = m); they don't work on submatrices yet. This *should* be a reasonably trivial thing to fix, as well.
//          -- Removed warning messages for generic calls
// 08.08.03 -- TRSM now works for all cases where SIDE == L and DIAG == N. DIAG == U is implemented but does not work correctly; SIDE == R is not yet implemented.
// 08.14.03 -- TRSM now works for all cases and accepts (and uses) leading-dimension information.
// 09.26.03 -- character input replaced with enumerated input to cause compiling errors and not run-time errors ( suggested by RAB ).

#ifndef _TEUCHOS_BLAS_HPP_
#define _TEUCHOS_BLAS_HPP_

/*! \file Teuchos_BLAS.hpp
    \brief Templated interface class to BLAS routines.
*/
/** \example BLAS/cxx_main.cpp
    This is an example of how to use the Teuchos::BLAS class.
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_BLAS_types.hpp"
#include "Teuchos_Assert.hpp"

/*! \class Teuchos::BLAS
    \brief Templated BLAS wrapper.

    The Teuchos::BLAS class provides functionality similar to the BLAS
    (Basic Linear Algebra Subprograms).  The BLAS provide portable,
    high- performance implementations of kernels such as dense vector
    sums, inner products, and norms (the BLAS 1 routines), dense
    matrix-vector multiplication and triangular solve (the BLAS 2
    routines), and dense matrix-matrix multiplication and triangular
    solve with multiple right-hand sides (the BLAS 3 routines).

    The standard BLAS interface is Fortran-specific.  Unfortunately,
    the interface between C++ and Fortran is not standard across all
    computer platforms.  The Teuchos::BLAS class provides C++ bindings
    for the BLAS kernels in order to insulate the rest of Petra from
    the details of C++ to Fortran translation.

    In addition to giving access to the standard BLAS functionality,
    Teuchos::BLAS also provides a generic fall-back implementation for
    any ScalarType class that defines the +, - * and / operators.

    Teuchos::BLAS only operates within a single shared-memory space,
    just like the BLAS.  It does not attempt to implement
    distributed-memory parallel matrix operations.

    \note This class has specializations for ScalarType=float and
      double, which use the BLAS library directly.  If you configure
      Teuchos to enable complex arithmetic support, via the CMake
      option -DTeuchos_ENABLE_COMPLEX:BOOL=ON, then this class will
      also invoke the BLAS library directly for
      ScalarType=std::complex<float> and std::complex<double>.
*/
namespace Teuchos
{
  extern TEUCHOSNUMERICS_LIB_DLL_EXPORT const char ESideChar[];
  extern TEUCHOSNUMERICS_LIB_DLL_EXPORT const char ETranspChar[];
  extern TEUCHOSNUMERICS_LIB_DLL_EXPORT const char EUploChar[];
  extern TEUCHOSNUMERICS_LIB_DLL_EXPORT const char EDiagChar[];
  extern TEUCHOSNUMERICS_LIB_DLL_EXPORT const char ETypeChar[];

  //! Default implementation for BLAS routines
  /*!
   * This class provides the default implementation for the BLAS routines.  It
   * is put in a separate class so that specializations of BLAS for other types
   * still have this implementation available.
   */
  template<typename OrdinalType, typename ScalarType>
  class DefaultBLASImpl
  {

    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;

  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    inline DefaultBLASImpl(void) {}

    //! Copy constructor.
    inline DefaultBLASImpl(const DefaultBLASImpl<OrdinalType, ScalarType>& /*BLAS_source*/) {}

    //! Destructor.
    inline virtual ~DefaultBLASImpl(void) {}
    //@}

    //! @name Level 1 BLAS Routines.
    //@{

    //! Computes a Givens plane rotation.
    void ROTG(ScalarType* da, ScalarType* db, MagnitudeType* c, ScalarType* s) const;

    //! Applies a Givens plane rotation.
    void ROT(const OrdinalType n, ScalarType* dx, const OrdinalType incx, ScalarType* dy, const OrdinalType incy, MagnitudeType* c, ScalarType* s) const;

    //! Scale the vector \c x by the constant \c alpha.
    void SCAL(const OrdinalType n, const ScalarType alpha, ScalarType* x, const OrdinalType incx) const;

    //! Copy the vector \c x to the vector \c y.
    void COPY(const OrdinalType n, const ScalarType* x, const OrdinalType incx, ScalarType* y, const OrdinalType incy) const;

    //! Perform the operation: \c y \c <- \c y+alpha*x.
    template <typename alpha_type, typename x_type>
    void AXPY(const OrdinalType n, const alpha_type alpha, const x_type* x, const OrdinalType incx, ScalarType* y, const OrdinalType incy) const;

    //! Sum the absolute values of the entries of \c x.
    typename ScalarTraits<ScalarType>::magnitudeType ASUM(const OrdinalType n, const ScalarType* x, const OrdinalType incx) const;

    //! Form the dot product of the vectors \c x and \c y.
    template <typename x_type, typename y_type>
    ScalarType DOT(const OrdinalType n, const x_type* x, const OrdinalType incx, const y_type* y, const OrdinalType incy) const;

    //! Compute the 2-norm of the vector \c x.
    typename ScalarTraits<ScalarType>::magnitudeType NRM2(const OrdinalType n, const ScalarType* x, const OrdinalType incx) const;

    //! Return the index of the element of \c x with the maximum magnitude.
    OrdinalType IAMAX(const OrdinalType n, const ScalarType* x, const OrdinalType incx) const;
    //@}

    //! @name Level 2 BLAS Routines.
    //@{

    //! Performs the matrix-vector operation:  \c y \c <- \c alpha*A*x+beta*y or \c y \c <- \c alpha*A'*x+beta*y where \c A is a general \c m by \c n matrix.
    template <typename alpha_type, typename A_type, typename x_type, typename beta_type>
    void GEMV(ETransp trans, const OrdinalType m, const OrdinalType n, const alpha_type alpha, const A_type* A,
              const OrdinalType lda, const x_type* x, const OrdinalType incx, const beta_type beta, ScalarType* y, const OrdinalType incy) const;

    //! Performs the matrix-vector operation:  \c x \c <- \c A*x or \c x \c <- \c A'*x where \c A is a unit/non-unit \c n by \c n upper/lower triangular matrix.
    template <typename A_type>
    void TRMV(EUplo uplo, ETransp trans, EDiag diag, const OrdinalType n, const A_type* A,
              const OrdinalType lda, ScalarType* x, const OrdinalType incx) const;

    //! \brief Performs the rank 1 operation:  \c A \c <- \c alpha*x*y'+A.
    /// \note  For complex arithmetic, this routine performs [Z/C]GERU.
    template <typename alpha_type, typename x_type, typename y_type>
    void GER(const OrdinalType m, const OrdinalType n, const alpha_type alpha, const x_type* x, const OrdinalType incx,
             const y_type* y, const OrdinalType incy, ScalarType* A, const OrdinalType lda) const;
    //@}

    //! @name Level 3 BLAS Routines.
    //@{

    /// \brief General matrix-matrix multiply.
    ///
    /// This computes C = alpha*op(A)*op(B) + beta*C.  op(X) here may
    /// be either X, the transpose of X, or the conjugate transpose of
    /// X.  op(A) has m rows and k columns, op(B) has k rows and n
    /// columns, and C has m rows and n columns.
    template <typename alpha_type, typename A_type, typename B_type, typename beta_type>
    void GEMM(ETransp transa, ETransp transb, const OrdinalType m, const OrdinalType n, const OrdinalType k, const alpha_type alpha, const A_type* A, const OrdinalType lda, const B_type* B, const OrdinalType ldb, const beta_type beta, ScalarType* C, const OrdinalType ldc) const;

    //! Performs the matrix-matrix operation: \c C \c <- \c alpha*A*B+beta*C or \c C \c <- \c alpha*B*A+beta*C where \c A is an \c m by \c m or \c n by \c n symmetric matrix and \c B is a general matrix.
    template <typename alpha_type, typename A_type, typename B_type, typename beta_type>
    void SYMM(ESide side, EUplo uplo, const OrdinalType m, const OrdinalType n, const alpha_type alpha, const A_type* A, const OrdinalType lda, const B_type* B, const OrdinalType ldb, const beta_type beta, ScalarType* C, const OrdinalType ldc) const;

    //! Performs the symmetric rank k operation: \c C <- \c alpha*A*A'+beta*C or \c C <- \c alpha*A'*A+beta*C, where \c alpha and \c beta are scalars, \c C is an \c n by \c n symmetric matrix and \c A is an \c n by \c k matrix in the first case or \c k by \c n matrix in the second case.
    template <typename alpha_type, typename A_type, typename beta_type>
    void SYRK(EUplo uplo, ETransp trans, const OrdinalType n, const OrdinalType k, const alpha_type alpha, const A_type* A, const OrdinalType lda, const beta_type beta, ScalarType* C, const OrdinalType ldc) const;

    //! Performs the matrix-matrix operation: \c B \c <- \c alpha*op(A)*B or \c B \c <- \c alpha*B*op(A) where \c op(A) is an unit/non-unit, upper/lower triangular matrix and \c B is a general matrix.
    template <typename alpha_type, typename A_type>
    void TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n,
                const alpha_type alpha, const A_type* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb) const;

    //! Solves the matrix equations:  \c op(A)*X=alpha*B or \c X*op(A)=alpha*B where \c X and \c B are \c m by \c n matrices, \c A is a unit/non-unit, upper/lower triangular matrix and \c op(A) is \c A or \c A'.  The matrix \c X is overwritten on \c B.
    template <typename alpha_type, typename A_type>
    void TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n,
                const alpha_type alpha, const A_type* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb) const;
    //@}
  };

  template<typename OrdinalType, typename ScalarType>
  class TEUCHOSNUMERICS_LIB_DLL_EXPORT BLAS : public DefaultBLASImpl<OrdinalType,ScalarType>
  {

    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;

  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    inline BLAS(void) {}

    //! Copy constructor.
    inline BLAS(const BLAS<OrdinalType, ScalarType>& /*BLAS_source*/) {}

    //! Destructor.
    inline virtual ~BLAS(void) {}
    //@}
  };

//------------------------------------------------------------------------------------------
//      LEVEL 1 BLAS ROUTINES
//------------------------------------------------------------------------------------------

  /// \namespace details
  /// \brief Teuchos implementation details.
  ///
  /// \warning Teuchos users should not use anything in this
  ///   namespace.  They should not even assume that the namespace
  ///   will continue to exist between releases.  The namespace's name
  ///   itself or anything it contains may change at any time.
  namespace details {

    // Compute magnitude.
    template<typename ScalarType, bool isComplex>
    class MagValue {
    public:
      void
      blas_dabs1(const ScalarType* a, typename ScalarTraits<ScalarType>::magnitudeType* ret) const; 
    };

    // Complex-arithmetic specialization. 
    template<typename ScalarType>
    class MagValue<ScalarType, true> {
    public:
      void
      blas_dabs1(const ScalarType* a, typename ScalarTraits<ScalarType>::magnitudeType* ret) const; 
    };

    // Real-arithmetic specialization.
    template<typename ScalarType>
    class MagValue<ScalarType, false> {
    public:
      void
      blas_dabs1(const ScalarType* a, ScalarType* ret) const;
    };
 
    template<typename ScalarType, bool isComplex>
    class GivensRotator {
    public:
      void
      ROTG (ScalarType* a,
            ScalarType* b,
            typename ScalarTraits<ScalarType>::magnitudeType* c,
            ScalarType* s) const;
    };

    // Complex-arithmetic specialization.
    template<typename ScalarType>
    class GivensRotator<ScalarType, true> {
    public:
      void
      ROTG (ScalarType* ca,
            ScalarType* cb,
            typename ScalarTraits<ScalarType>::magnitudeType* c,
            ScalarType* s) const;
    };

    // Real-arithmetic specialization.
    template<typename ScalarType>
    class GivensRotator<ScalarType, false> {
    public:
      void
      ROTG (ScalarType* da,
            ScalarType* db,
            ScalarType* c,
            ScalarType* s) const;
    private:
      /// Return ABS(x) if y > 0 or y is +0, else -ABS(x) (if y is -0 or < 0).
      ///
      /// Note that SIGN respects IEEE 754 floating-point signed zero.
      /// This is a hopefully correct implementation of the Fortran
      /// type-generic SIGN intrinsic.  ROTG for complex arithmetic
      /// doesn't require this function.  C99 provides a copysign()
      /// math library function, but we are not able to rely on the
      /// existence of C99 functions here.
      ///
      /// We provide this method on purpose only for the
      /// real-arithmetic specialization of GivensRotator.  Complex
      /// numbers don't have a sign; they have an angle.
      ScalarType SIGN (ScalarType x, ScalarType y) const {
        typedef ScalarTraits<ScalarType> STS;

        if (y > STS::zero()) {
          return STS::magnitude (x);
        } else if (y < STS::zero()) {
          return -STS::magnitude (x);
        } else { // y == STS::zero()
          // Suppose that ScalarType implements signed zero, as IEEE
          // 754 - compliant floating-point numbers should.  You can't
          // use == to test for signed zero, since +0 == -0.  However,
          // 1/0 = Inf > 0 and 1/-0 = -Inf < 0.  Let's hope ScalarType
          // supports Inf... we don't need to test for Inf, just see
          // if it's greater than or less than zero.
          //
          // NOTE: This ONLY works if ScalarType is real.  Complex
          // infinity doesn't have a sign, so we can't compare it with
          // zero.  That's OK, because finite complex numbers don't
          // have a sign either; they have an angle.
          ScalarType signedInfinity = STS::one() / y;
          if (signedInfinity > STS::zero()) {
            return STS::magnitude (x);
          } else {
            // Even if ScalarType doesn't implement signed zero,
            // Fortran's SIGN intrinsic returns -ABS(X) if the second
            // argument Y is zero.  We imitate this behavior here.
            return -STS::magnitude (x);
          }
        }
      }
    };

    // Implementation of complex-arithmetic specialization.
    template<typename ScalarType>
    void
    GivensRotator<ScalarType, true>::
    ROTG (ScalarType* ca,
          ScalarType* cb,
          typename ScalarTraits<ScalarType>::magnitudeType* c,
          ScalarType* s) const
    {
      typedef ScalarTraits<ScalarType> STS;
      typedef typename STS::magnitudeType MagnitudeType;
      typedef ScalarTraits<MagnitudeType> STM;

      // This is a straightforward translation into C++ of the
      // reference BLAS' implementation of ZROTG.  You can get
      // the Fortran 77 source code of ZROTG here:
      //
      // http://www.netlib.org/blas/zrotg.f
      //
      // I used the following rules to translate Fortran types and
      // intrinsic functions into C++:
      //
      // DOUBLE PRECISION -> MagnitudeType
      // DOUBLE COMPLEX -> ScalarType
      // CDABS -> STS::magnitude
      // DCMPLX -> ScalarType constructor (assuming that ScalarType
      //   is std::complex<MagnitudeType>)
      // DCONJG -> STS::conjugate
      // DSQRT -> STM::squareroot
      ScalarType alpha;
      MagnitudeType norm, scale;

      if (STS::magnitude (*ca) == STM::zero()) {
        *c = STM::zero();
        *s = STS::one();
        *ca = *cb;
      } else {
        scale = STS::magnitude (*ca) + STS::magnitude (*cb);
        { // I introduced temporaries into the translated BLAS code in
          // order to make the expression easier to read and also save a
          // few floating-point operations.
          const MagnitudeType ca_scaled =
            STS::magnitude (*ca / ScalarType(scale, STM::zero()));
          const MagnitudeType cb_scaled =
            STS::magnitude (*cb / ScalarType(scale, STM::zero()));
          norm = scale *
            STM::squareroot (ca_scaled*ca_scaled + cb_scaled*cb_scaled);
        }
        alpha = *ca / STS::magnitude (*ca);
        *c = STS::magnitude (*ca) / norm;
        *s = alpha * STS::conjugate (*cb) / norm;
        *ca = alpha * norm;
      }
    }

    // Implementation of real-arithmetic specialization.
    template<typename ScalarType>
    void
    GivensRotator<ScalarType, false>::
    ROTG (ScalarType* da,
          ScalarType* db,
          ScalarType* c,
          ScalarType* s) const
    {
      typedef ScalarTraits<ScalarType> STS;

      // This is a straightforward translation into C++ of the
      // reference BLAS' implementation of DROTG.  You can get
      // the Fortran 77 source code of DROTG here:
      //
      // http://www.netlib.org/blas/drotg.f
      //
      // I used the following rules to translate Fortran types and
      // intrinsic functions into C++:
      //
      // DOUBLE PRECISION -> ScalarType
      // DABS -> STS::magnitude
      // DSQRT -> STM::squareroot
      // DSIGN -> SIGN (see below)
      //
      // DSIGN(x,y) (the old DOUBLE PRECISION type-specific form of
      // the Fortran type-generic SIGN intrinsic) required special
      // translation, which we did in a separate utility function in
      // the specializaton of GivensRotator for real arithmetic.
      // (ROTG for complex arithmetic doesn't require this function.)
      // C99 provides a copysign() math library function, but we are
      // not able to rely on the existence of C99 functions here.
      ScalarType r, roe, scale, z;

      roe = *db;
      if (STS::magnitude (*da) > STS::magnitude (*db)) {
        roe = *da;
      }
      scale = STS::magnitude (*da) + STS::magnitude (*db);
      if (scale == STS::zero()) {
        *c = STS::one();
        *s = STS::zero();
        r = STS::zero();
        z = STS::zero();
      } else {
        // I introduced temporaries into the translated BLAS code in
        // order to make the expression easier to read and also save
        // a few floating-point operations.
        const ScalarType da_scaled = *da / scale;
        const ScalarType db_scaled = *db / scale;
        r = scale * STS::squareroot (da_scaled*da_scaled + db_scaled*db_scaled);
        r = SIGN (STS::one(), roe) * r;
        *c = *da / r;
        *s = *db / r;
        z = STS::one();
        if (STS::magnitude (*da) > STS::magnitude (*db)) {
          z = *s;
        }
        if (STS::magnitude (*db) >= STS::magnitude (*da) && *c != STS::zero()) {
          z = STS::one() / *c;
        }
      }

      *da = r;
      *db = z;
    }

    // Real-valued implementation of MagValue
    template<typename ScalarType>
    void 
    MagValue<ScalarType, false>::
    blas_dabs1(const ScalarType* a, ScalarType* ret) const 
    {
      *ret = Teuchos::ScalarTraits<ScalarType>::magnitude( *a );  
    }

    // Complex-valued implementation of MagValue
    template<typename ScalarType>
    void 
    MagValue<ScalarType, true>::
    blas_dabs1(const ScalarType* a, typename ScalarTraits<ScalarType>::magnitudeType* ret) const 
    {
      *ret = ScalarTraits<typename ScalarTraits<ScalarType>::magnitudeType>::magnitude(a->real());
      *ret += ScalarTraits<typename ScalarTraits<ScalarType>::magnitudeType>::magnitude(a->imag());
    }

  } // namespace details

  template<typename OrdinalType, typename ScalarType>
  void
  DefaultBLASImpl<OrdinalType, ScalarType>::
  ROTG (ScalarType* da,
        ScalarType* db,
        MagnitudeType* c,
        ScalarType* s) const
  {
    typedef ScalarTraits<ScalarType> STS;
    details::GivensRotator<ScalarType, STS::isComplex> rotator;
    rotator.ROTG (da, db, c, s);
  }

  template<typename OrdinalType, typename ScalarType>
  void DefaultBLASImpl<OrdinalType,ScalarType>::ROT(const OrdinalType n, ScalarType* dx, const OrdinalType incx, ScalarType* dy, const OrdinalType incy, MagnitudeType* c, ScalarType* s) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    ScalarType sconj = Teuchos::ScalarTraits<ScalarType>::conjugate(*s);
    if (n <= 0) return;
    if (incx==1 && incy==1) {
      for(OrdinalType i=0; i<n; ++i) {
        ScalarType temp = *c*dx[i] + sconj*dy[i];
        dy[i] = *c*dy[i] - sconj*dx[i];
        dx[i] = temp;
      }
    }
    else {
      OrdinalType ix = 0, iy = 0;
      if (incx < izero) ix = (-n+1)*incx;
      if (incy < izero) iy = (-n+1)*incy;
      for(OrdinalType i=0; i<n; ++i) {
        ScalarType temp = *c*dx[ix] + sconj*dy[iy];
        dy[iy] = *c*dy[iy] - sconj*dx[ix];
        dx[ix] = temp;
        ix += incx;
        iy += incy;
      }
    }
  }

  template<typename OrdinalType, typename ScalarType>
  void DefaultBLASImpl<OrdinalType, ScalarType>::SCAL(const OrdinalType n, const ScalarType alpha, ScalarType* x, const OrdinalType incx) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    OrdinalType i, ix = izero;

    if ( n < ione || incx < ione )
      return;

    // Scale the vector.
    for(i = izero; i < n; i++)
    {
      x[ix] *= alpha;
      ix += incx;
    }
  } /* end SCAL */

  template<typename OrdinalType, typename ScalarType>
  void DefaultBLASImpl<OrdinalType, ScalarType>::COPY(const OrdinalType n, const ScalarType* x, const OrdinalType incx, ScalarType* y, const OrdinalType incy) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    OrdinalType i, ix = izero, iy = izero;
    if ( n > izero ) {
        // Set the initial indices (ix, iy).
        if (incx < izero) { ix = (-n+ione)*incx; }
        if (incy < izero) { iy = (-n+ione)*incy; }

        for(i = izero; i < n; i++)
          {
            y[iy] = x[ix];
            ix += incx;
            iy += incy;
          }
      }
  } /* end COPY */

  template<typename OrdinalType, typename ScalarType>
  template <typename alpha_type, typename x_type>
  void DefaultBLASImpl<OrdinalType, ScalarType>::AXPY(const OrdinalType n, const alpha_type alpha, const x_type* x, const OrdinalType incx, ScalarType* y, const OrdinalType incy) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    OrdinalType i, ix = izero, iy = izero;
    if( n > izero && alpha != ScalarTraits<alpha_type>::zero())
      {
        // Set the initial indices (ix, iy).
        if (incx < izero) { ix = (-n+ione)*incx; }
        if (incy < izero) { iy = (-n+ione)*incy; }

        for(i = izero; i < n; i++)
          {
            y[iy] += alpha * x[ix];
            ix += incx;
            iy += incy;
          }
      }
  } /* end AXPY */

  template<typename OrdinalType, typename ScalarType>
  typename ScalarTraits<ScalarType>::magnitudeType DefaultBLASImpl<OrdinalType, ScalarType>::ASUM(const OrdinalType n, const ScalarType* x, const OrdinalType incx) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    typename ScalarTraits<ScalarType>::magnitudeType temp, result =
      ScalarTraits<typename ScalarTraits<ScalarType>::magnitudeType>::zero();
    OrdinalType i, ix = izero;

    if ( n < ione || incx < ione )
      return result;

    details::MagValue<ScalarType, ScalarTraits<ScalarType>::isComplex> mval;
    for (i = izero; i < n; i++)
      {
        mval.blas_dabs1( &x[ix], &temp );
        result += temp;
        ix += incx;
      }
   
    return result;
  } /* end ASUM */

  template<typename OrdinalType, typename ScalarType>
  template <typename x_type, typename y_type>
  ScalarType DefaultBLASImpl<OrdinalType, ScalarType>::DOT(const OrdinalType n, const x_type* x, const OrdinalType incx, const y_type* y, const OrdinalType incy) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType result = ScalarTraits<ScalarType>::zero();
    OrdinalType i, ix = izero, iy = izero;
    if( n > izero )
      {
        // Set the initial indices (ix,iy).
        if (incx < izero) { ix = (-n+ione)*incx; }
        if (incy < izero) { iy = (-n+ione)*incy; }

        for(i = izero; i < n; i++)
          {
            result += ScalarTraits<x_type>::conjugate(x[ix]) * y[iy];
            ix += incx;
            iy += incy;
          }
      }
    return result;
  } /* end DOT */

  template<typename OrdinalType, typename ScalarType>
  typename ScalarTraits<ScalarType>::magnitudeType DefaultBLASImpl<OrdinalType, ScalarType>::NRM2(const OrdinalType n, const ScalarType* x, const OrdinalType incx) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    typename ScalarTraits<ScalarType>::magnitudeType result =
      ScalarTraits<typename ScalarTraits<ScalarType>::magnitudeType>::zero();
    OrdinalType i, ix = izero;

    if ( n < ione || incx < ione )
      return result;

    for(i = izero; i < n; i++)
      {
        result += ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::conjugate(x[ix]) * x[ix]);
        ix += incx;
      }
    result = ScalarTraits<typename ScalarTraits<ScalarType>::magnitudeType>::squareroot(result);
    return result;
  } /* end NRM2 */

  template<typename OrdinalType, typename ScalarType>
  OrdinalType DefaultBLASImpl<OrdinalType, ScalarType>::IAMAX(const OrdinalType n, const ScalarType* x, const OrdinalType incx) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    OrdinalType result = izero, ix = izero, i;
    typename ScalarTraits<ScalarType>::magnitudeType curval =
      ScalarTraits<typename ScalarTraits<ScalarType>::magnitudeType>::zero();
    typename ScalarTraits<ScalarType>::magnitudeType maxval =
      ScalarTraits<typename ScalarTraits<ScalarType>::magnitudeType>::zero();

    if ( n < ione || incx < ione )
      return result;
  
    details::MagValue<ScalarType, ScalarTraits<ScalarType>::isComplex> mval;

    mval.blas_dabs1( &x[ix], &maxval );
    ix += incx;
    for(i = ione; i < n; i++)
      {
        mval.blas_dabs1( &x[ix], &curval );
        if(curval > maxval)
          {
            result = i;
            maxval = curval;
          }
        ix += incx;
      }

    return result + 1; // the BLAS I?AMAX functions return 1-indexed (Fortran-style) values
  } /* end IAMAX */

//------------------------------------------------------------------------------------------
//      LEVEL 2 BLAS ROUTINES
//------------------------------------------------------------------------------------------
  template<typename OrdinalType, typename ScalarType>
  template <typename alpha_type, typename A_type, typename x_type, typename beta_type>
  void DefaultBLASImpl<OrdinalType, ScalarType>::GEMV(ETransp trans, const OrdinalType m, const OrdinalType n, const alpha_type alpha, const A_type* A, const OrdinalType lda, const x_type* x, const OrdinalType incx, const beta_type beta, ScalarType* y, const OrdinalType incy) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    alpha_type alpha_zero = ScalarTraits<alpha_type>::zero();
    beta_type beta_zero = ScalarTraits<beta_type>::zero();
    x_type x_zero = ScalarTraits<x_type>::zero();
    ScalarType y_zero = ScalarTraits<ScalarType>::zero();
    beta_type beta_one = ScalarTraits<beta_type>::one();
    bool noConj = true;
    bool BadArgument = false;

    // Quick return if there is nothing to do!
    if( m == izero || n == izero || ( alpha == alpha_zero && beta == beta_one ) ){ return; }

    // Otherwise, we need to check the argument list.
    if( m < izero ) {
        std::cout << "BLAS::GEMV Error: M == " << m << std::endl;
        BadArgument = true;
    }
    if( n < izero ) {
        std::cout << "BLAS::GEMV Error: N == " << n << std::endl;
        BadArgument = true;
    }
    if( lda < m ) {
        std::cout << "BLAS::GEMV Error: LDA < MAX(1,M)"<< std::endl;
        BadArgument = true;
    }
    if( incx == izero ) {
        std::cout << "BLAS::GEMV Error: INCX == 0"<< std::endl;
        BadArgument = true;
    }
    if( incy == izero ) {
        std::cout << "BLAS::GEMV Error: INCY == 0"<< std::endl;
        BadArgument = true;
    }

    if(!BadArgument) {
      OrdinalType i, j, lenx, leny, ix, iy, jx, jy;
      OrdinalType kx = izero, ky = izero;
      ScalarType temp;

      // Determine the lengths of the vectors x and y.
      if(ETranspChar[trans] == 'N') {
        lenx = n;
        leny = m;
      } else {
        lenx = m;
        leny = n;
      }

      // Determine if this is a conjugate tranpose
      noConj = (ETranspChar[trans] == 'T');

      // Set the starting pointers for the vectors x and y if incx/y < 0.
      if (incx < izero ) { kx =  (ione - lenx)*incx; }
      if (incy < izero ) { ky =  (ione - leny)*incy; }

      // Form y = beta*y
      ix = kx; iy = ky;
      if(beta != beta_one) {
        if (incy == ione) {
          if (beta == beta_zero) {
            for(i = izero; i < leny; i++) { y[i] = y_zero; }
          } else {
            for(i = izero; i < leny; i++) { y[i] *= beta; }
          }
        } else {
          if (beta == beta_zero) {
            for(i = izero; i < leny; i++) {
              y[iy] = y_zero;
              iy += incy;
            }
          } else {
            for(i = izero; i < leny; i++) {
              y[iy] *= beta;
              iy += incy;
            }
          }
        }
      }

      // Return if we don't have to do anything more.
      if(alpha == alpha_zero) { return; }

      if( ETranspChar[trans] == 'N' ) {
        // Form y = alpha*A*y
        jx = kx;
        if (incy == ione) {
          for(j = izero; j < n; j++) {
            if (x[jx] != x_zero) {
              temp = alpha*x[jx];
              for(i = izero; i < m; i++) {
                y[i] += temp*A[j*lda + i];
              }
            }
            jx += incx;
          }
        } else {
          for(j = izero; j < n; j++) {
            if (x[jx] != x_zero) {
              temp = alpha*x[jx];
              iy = ky;
              for(i = izero; i < m; i++) {
                y[iy] += temp*A[j*lda + i];
                iy += incy;
              }
            }
            jx += incx;
          }
        }
      } else {
        jy = ky;
        if (incx == ione) {
          for(j = izero; j < n; j++) {
            temp = y_zero;
            if ( noConj ) {
              for(i = izero; i < m; i++) {
                temp += A[j*lda + i]*x[i];
              }
            } else {
              for(i = izero; i < m; i++) {
                temp += ScalarTraits<A_type>::conjugate(A[j*lda + i])*x[i];
              }
            }
            y[jy] += alpha*temp;
            jy += incy;
          }
        } else {
          for(j = izero; j < n; j++) {
            temp = y_zero;
            ix = kx;
            if ( noConj ) {
              for (i = izero; i < m; i++) {
                temp += A[j*lda + i]*x[ix];
                ix += incx;
              }
            } else {
              for (i = izero; i < m; i++) {
                temp += ScalarTraits<A_type>::conjugate(A[j*lda + i])*x[ix];
                ix += incx;
              }
            }
            y[jy] += alpha*temp;
            jy += incy;
          }
        }
      }
    } /* if (!BadArgument) */
  } /* end GEMV */

 template<typename OrdinalType, typename ScalarType>
 template <typename A_type>
 void DefaultBLASImpl<OrdinalType, ScalarType>::TRMV(EUplo uplo, ETransp trans, EDiag diag, const OrdinalType n, const A_type* A, const OrdinalType lda, ScalarType* x, const OrdinalType incx) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    ScalarType zero = ScalarTraits<ScalarType>::zero();
    bool BadArgument = false;
    bool noConj = true;

    // Quick return if there is nothing to do!
    if( n == izero ){ return; }

    // Otherwise, we need to check the argument list.
    if( n < izero ) {
      std::cout << "BLAS::TRMV Error: N == " << n << std::endl;
      BadArgument = true;
    }
    if( lda < n ) {
      std::cout << "BLAS::TRMV Error: LDA < MAX(1,N)"<< std::endl;
      BadArgument = true;
    }
    if( incx == izero ) {
      std::cout << "BLAS::TRMV Error: INCX == 0"<< std::endl;
      BadArgument = true;
    }

    if(!BadArgument) {
      OrdinalType i, j, ix, jx, kx = izero;
      ScalarType temp;
      bool noUnit = (EDiagChar[diag] == 'N');

      // Determine if this is a conjugate tranpose
      noConj = (ETranspChar[trans] == 'T');

      // Set the starting pointer for the vector x if incx < 0.
      if (incx < izero) { kx = (-n+ione)*incx; }

      // Start the operations for a nontransposed triangular matrix
      if (ETranspChar[trans] == 'N') {
        /* Compute x = A*x */
        if (EUploChar[uplo] == 'U') {
          /* A is an upper triangular matrix */
          if (incx == ione) {
            for (j=izero; j<n; j++) {
              if (x[j] != zero) {
                temp = x[j];
                for (i=izero; i<j; i++) {
                  x[i] += temp*A[j*lda + i];
                }
                if ( noUnit )
                  x[j] *= A[j*lda + j];
              }
            }
          } else {
            jx = kx;
            for (j=izero; j<n; j++) {
              if (x[jx] != zero) {
                temp = x[jx];
                ix = kx;
                for (i=izero; i<j; i++) {
                  x[ix] += temp*A[j*lda + i];
                  ix += incx;
                }
                if ( noUnit )
                  x[jx] *= A[j*lda + j];
              }
              jx += incx;
            }
          } /* if (incx == ione) */
        } else { /* A is a lower triangular matrix */
          if (incx == ione) {
            for (j=n-ione; j>-ione; j--) {
              if (x[j] != zero) {
                temp = x[j];
                for (i=n-ione; i>j; i--) {
                  x[i] += temp*A[j*lda + i];
                }
                if ( noUnit )
                  x[j] *= A[j*lda + j];
              }
            }
          } else {
            kx += (n-ione)*incx;
            jx = kx;
            for (j=n-ione; j>-ione; j--) {
              if (x[jx] != zero) {
                temp = x[jx];
                ix = kx;
                for (i=n-ione; i>j; i--) {
                  x[ix] += temp*A[j*lda + i];
                  ix -= incx;
                }
                if ( noUnit )
                  x[jx] *= A[j*lda + j];
              }
              jx -= incx;
            }
          }
        } /* if (EUploChar[uplo]=='U') */
      } else { /* A is transposed/conjugated */
        /* Compute x = A'*x */
        if (EUploChar[uplo]=='U') {
          /* A is an upper triangular matrix */
          if (incx == ione) {
            for (j=n-ione; j>-ione; j--) {
              temp = x[j];
              if ( noConj ) {
                if ( noUnit )
                  temp *= A[j*lda + j];
                for (i=j-ione; i>-ione; i--) {
                  temp += A[j*lda + i]*x[i];
                }
              } else {
                if ( noUnit )
                  temp *= ScalarTraits<A_type>::conjugate(A[j*lda + j]);
                for (i=j-ione; i>-ione; i--) {
                  temp += ScalarTraits<A_type>::conjugate(A[j*lda + i])*x[i];
                }
              }
              x[j] = temp;
            }
          } else {
            jx = kx + (n-ione)*incx;
            for (j=n-ione; j>-ione; j--) {
              temp = x[jx];
              ix = jx;
              if ( noConj ) {
                if ( noUnit )
                  temp *= A[j*lda + j];
                for (i=j-ione; i>-ione; i--) {
                  ix -= incx;
                  temp += A[j*lda + i]*x[ix];
                }
              } else {
                if ( noUnit )
                  temp *= ScalarTraits<A_type>::conjugate(A[j*lda + j]);
                for (i=j-ione; i>-ione; i--) {
                  ix -= incx;
                  temp += ScalarTraits<A_type>::conjugate(A[j*lda + i])*x[ix];
                }
              }
              x[jx] = temp;
              jx -= incx;
            }
          }
        } else {
          /* A is a lower triangular matrix */
          if (incx == ione) {
            for (j=izero; j<n; j++) {
              temp = x[j];
              if ( noConj ) {
                if ( noUnit )
                  temp *= A[j*lda + j];
                for (i=j+ione; i<n; i++) {
                  temp += A[j*lda + i]*x[i];
                }
              } else {
                if ( noUnit )
                  temp *= ScalarTraits<A_type>::conjugate(A[j*lda + j]);
                for (i=j+ione; i<n; i++) {
                  temp += ScalarTraits<A_type>::conjugate(A[j*lda + i])*x[i];
                }
              }
              x[j] = temp;
            }
          } else {
            jx = kx;
            for (j=izero; j<n; j++) {
              temp = x[jx];
              ix = jx;
              if ( noConj ) {
                if ( noUnit )
                  temp *= A[j*lda + j];
                for (i=j+ione; i<n; i++) {
                  ix += incx;
                  temp += A[j*lda + i]*x[ix];
                }
              } else {
                if ( noUnit )
                  temp *= ScalarTraits<A_type>::conjugate(A[j*lda + j]);
                for (i=j+ione; i<n; i++) {
                  ix += incx;
                  temp += ScalarTraits<A_type>::conjugate(A[j*lda + i])*x[ix];
                }
              }
              x[jx] = temp;
              jx += incx;
            }
          }
        } /* if (EUploChar[uplo]=='U') */
      } /* if (ETranspChar[trans]=='N') */
    } /* if (!BadArgument) */
  } /* end TRMV */

  template<typename OrdinalType, typename ScalarType>
  template <typename alpha_type, typename x_type, typename y_type>
  void DefaultBLASImpl<OrdinalType, ScalarType>::GER(const OrdinalType m, const OrdinalType n, const alpha_type alpha, const x_type* x, const OrdinalType incx, const y_type* y, const OrdinalType incy, ScalarType* A, const OrdinalType lda) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    alpha_type alpha_zero = ScalarTraits<alpha_type>::zero();
    y_type y_zero = ScalarTraits<y_type>::zero();
    bool BadArgument = false;

    // Quick return if there is nothing to do!
    if( m == izero || n == izero || alpha == alpha_zero ){ return; }

    // Otherwise, we need to check the argument list.
    if( m < izero ) {
        std::cout << "BLAS::GER Error: M == " << m << std::endl;
        BadArgument = true;
    }
    if( n < izero ) {
        std::cout << "BLAS::GER Error: N == " << n << std::endl;
        BadArgument = true;
    }
    if( lda < m ) {
        std::cout << "BLAS::GER Error: LDA < MAX(1,M)"<< std::endl;
        BadArgument = true;
    }
    if( incx == 0 ) {
        std::cout << "BLAS::GER Error: INCX == 0"<< std::endl;
        BadArgument = true;
    }
    if( incy == 0 ) {
        std::cout << "BLAS::GER Error: INCY == 0"<< std::endl;
        BadArgument = true;
    }

    if(!BadArgument) {
      OrdinalType i, j, ix, jy = izero, kx = izero;
      ScalarType temp;

      // Set the starting pointers for the vectors x and y if incx/y < 0.
      if (incx < izero) { kx = (-m+ione)*incx; }
      if (incy < izero) { jy = (-n+ione)*incy; }

      // Start the operations for incx == 1
      if( incx == ione ) {
        for( j=izero; j<n; j++ ) {
          if ( y[jy] != y_zero ) {
            temp = alpha*y[jy];
            for ( i=izero; i<m; i++ ) {
              A[j*lda + i] += x[i]*temp;
            }
          }
          jy += incy;
        }
      }
      else { // Start the operations for incx != 1
        for( j=izero; j<n; j++ ) {
          if ( y[jy] != y_zero ) {
            temp = alpha*y[jy];
            ix = kx;
            for( i=izero; i<m; i++ ) {
              A[j*lda + i] += x[ix]*temp;
              ix += incx;
            }
          }
          jy += incy;
        }
      }
    } /* if(!BadArgument) */
  } /* end GER */

//------------------------------------------------------------------------------------------
//      LEVEL 3 BLAS ROUTINES
//------------------------------------------------------------------------------------------

  template<typename OrdinalType, typename ScalarType>
  template <typename alpha_type, typename A_type, typename B_type, typename beta_type>
  void DefaultBLASImpl<OrdinalType, ScalarType>::GEMM(ETransp transa, ETransp transb, const OrdinalType m, const OrdinalType n, const OrdinalType k, const alpha_type alpha, const A_type* A, const OrdinalType lda, const B_type* B, const OrdinalType ldb, const beta_type beta, ScalarType* C, const OrdinalType ldc) const
  {

    typedef TypeNameTraits<OrdinalType> OTNT;
    typedef TypeNameTraits<ScalarType> STNT;

    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    alpha_type alpha_zero = ScalarTraits<alpha_type>::zero();
    beta_type beta_zero = ScalarTraits<beta_type>::zero();
    B_type B_zero = ScalarTraits<B_type>::zero();
    ScalarType C_zero = ScalarTraits<ScalarType>::zero();
    beta_type beta_one = ScalarTraits<beta_type>::one();
    OrdinalType i, j, p;
    OrdinalType NRowA = m, NRowB = k;
    ScalarType temp;
    bool BadArgument = false;
    bool conjA = false, conjB = false;

    // Change dimensions of matrix if either matrix is transposed
    if( !(ETranspChar[transa]=='N') ) {
      NRowA = k;
    }
    if( !(ETranspChar[transb]=='N') ) {
      NRowB = n;
    }

    // Quick return if there is nothing to do!
    if( (m==izero) || (n==izero) || (((alpha==alpha_zero)||(k==izero)) && (beta==beta_one)) ){ return; }
    if( m < izero ) {
      std::cout << "BLAS::GEMM Error: M == " << m << std::endl;
      BadArgument = true;
    }
    if( n < izero ) {
      std::cout << "BLAS::GEMM Error: N == " << n << std::endl;
      BadArgument = true;
    }
    if( k < izero ) {
      std::cout << "BLAS::GEMM Error: K == " << k << std::endl;
      BadArgument = true;
    }
    if( lda < NRowA ) {
      std::cout << "BLAS::GEMM Error: LDA < "<<NRowA<<std::endl;
      BadArgument = true;
    }
    if( ldb < NRowB ) {
      std::cout << "BLAS::GEMM Error: LDB < "<<NRowB<<std::endl;
      BadArgument = true;
    }
     if( ldc < m ) {
      std::cout << "BLAS::GEMM Error: LDC < MAX(1,M)"<< std::endl;
      BadArgument = true;
    }

    if(!BadArgument) {

      // Determine if this is a conjugate tranpose
      conjA = (ETranspChar[transa] == 'C');
      conjB = (ETranspChar[transb] == 'C');

      // Only need to scale the resulting matrix C.
      if( alpha == alpha_zero ) {
        if( beta == beta_zero ) {
          for (j=izero; j<n; j++) {
            for (i=izero; i<m; i++) {
              C[j*ldc + i] = C_zero;
            }
          }
        } else {
          for (j=izero; j<n; j++) {
            for (i=izero; i<m; i++) {
              C[j*ldc + i] *= beta;
            }
          }
        }
        return;
      }
      //
      // Now start the operations.
      //
      if ( ETranspChar[transb]=='N' ) {
        if ( ETranspChar[transa]=='N' ) {
          // Compute C = alpha*A*B + beta*C
          for (j=izero; j<n; j++) {
            if( beta == beta_zero ) {
              for (i=izero; i<m; i++){
                C[j*ldc + i] = C_zero;
              }
            } else if( beta != beta_one ) {
              for (i=izero; i<m; i++){
                C[j*ldc + i] *= beta;
              }
            }
            for (p=izero; p<k; p++){
              if (B[j*ldb + p] != B_zero ){
                temp = alpha*B[j*ldb + p];
                for (i=izero; i<m; i++) {
                  C[j*ldc + i] += temp*A[p*lda + i];
                }
              }
            }
          }
        } else if ( conjA ) {
          // Compute C = alpha*conj(A')*B + beta*C
          for (j=izero; j<n; j++) {
            for (i=izero; i<m; i++) {
              temp = C_zero;
              for (p=izero; p<k; p++) {
                temp += ScalarTraits<A_type>::conjugate(A[i*lda + p])*B[j*ldb + p];
              }
              if (beta == beta_zero) {
                C[j*ldc + i] = alpha*temp;
              } else {
                C[j*ldc + i] = alpha*temp + beta*C[j*ldc + i];
              }
            }
          }
        } else {
          // Compute C = alpha*A'*B + beta*C
          for (j=izero; j<n; j++) {
            for (i=izero; i<m; i++) {
              temp = C_zero;
              for (p=izero; p<k; p++) {
                temp += A[i*lda + p]*B[j*ldb + p];
              }
              if (beta == beta_zero) {
                C[j*ldc + i] = alpha*temp;
              } else {
                C[j*ldc + i] = alpha*temp + beta*C[j*ldc + i];
              }
            }
          }
        }
      } else if ( ETranspChar[transa]=='N' ) {
        if ( conjB ) {
          // Compute C = alpha*A*conj(B') + beta*C
          for (j=izero; j<n; j++) {
            if (beta == beta_zero) {
              for (i=izero; i<m; i++) {
                C[j*ldc + i] = C_zero;
              }
            } else if ( beta != beta_one ) {
              for (i=izero; i<m; i++) {
                C[j*ldc + i] *= beta;
              }
            }
            for (p=izero; p<k; p++) {
              if (B[p*ldb + j] != B_zero) {
                temp = alpha*ScalarTraits<B_type>::conjugate(B[p*ldb + j]);
                for (i=izero; i<m; i++) {
                  C[j*ldc + i] += temp*A[p*lda + i];
                }
              }
            }
          }
        } else {
          // Compute C = alpha*A*B' + beta*C
          for (j=izero; j<n; j++) {
            if (beta == beta_zero) {
              for (i=izero; i<m; i++) {
                C[j*ldc + i] = C_zero;
              }
            } else if ( beta != beta_one ) {
              for (i=izero; i<m; i++) {
                C[j*ldc + i] *= beta;
              }
            }
            for (p=izero; p<k; p++) {
              if (B[p*ldb + j] != B_zero) {
                temp = alpha*B[p*ldb + j];
                for (i=izero; i<m; i++) {
                  C[j*ldc + i] += temp*A[p*lda + i];
                }
              }
            }
          }
        }
      } else if ( conjA ) {
        if ( conjB ) {
          // Compute C = alpha*conj(A')*conj(B') + beta*C
          for (j=izero; j<n; j++) {
            for (i=izero; i<m; i++) {
              temp = C_zero;
              for (p=izero; p<k; p++) {
                temp += ScalarTraits<A_type>::conjugate(A[i*lda + p])
                      * ScalarTraits<B_type>::conjugate(B[p*ldb + j]);
              }
              if (beta == beta_zero) {
                C[j*ldc + i] = alpha*temp;
              } else {
                C[j*ldc + i] = alpha*temp + beta*C[j*ldc + i];
              }
            }
          }
        } else {
          // Compute C = alpha*conj(A')*B' + beta*C
          for (j=izero; j<n; j++) {
            for (i=izero; i<m; i++) {
              temp = C_zero;
              for (p=izero; p<k; p++) {
                temp += ScalarTraits<A_type>::conjugate(A[i*lda + p])
                      * B[p*ldb + j];
              }
              if (beta == beta_zero) {
                C[j*ldc + i] = alpha*temp;
              } else {
                C[j*ldc + i] = alpha*temp + beta*C[j*ldc + i];
              }
            }
          }
        }
     } else {
       if ( conjB ) {
         // Compute C = alpha*A'*conj(B') + beta*C
         for (j=izero; j<n; j++) {
            for (i=izero; i<m; i++) {
              temp = C_zero;
              for (p=izero; p<k; p++) {
                temp += A[i*lda + p]
                      * ScalarTraits<B_type>::conjugate(B[p*ldb + j]);
              }
              if (beta == beta_zero) {
                C[j*ldc + i] = alpha*temp;
              } else {
                C[j*ldc + i] = alpha*temp + beta*C[j*ldc + i];
              }
            }
          }
        } else {
          // Compute C = alpha*A'*B' + beta*C
          for (j=izero; j<n; j++) {
            for (i=izero; i<m; i++) {
              temp = C_zero;
              for (p=izero; p<k; p++) {
                temp += A[i*lda + p]*B[p*ldb + j];
              }
              if (beta == beta_zero) {
                C[j*ldc + i] = alpha*temp;
              } else {
                C[j*ldc + i] = alpha*temp + beta*C[j*ldc + i];
              }
            }
          }
        } // end if (ETranspChar[transa]=='N') ...
      } // end if (ETranspChar[transb]=='N') ...
    } // end if (!BadArgument) ...
  } // end of GEMM


  template<typename OrdinalType, typename ScalarType>
  template <typename alpha_type, typename A_type, typename B_type, typename beta_type>
  void DefaultBLASImpl<OrdinalType, ScalarType>::SYMM(ESide side, EUplo uplo, const OrdinalType m, const OrdinalType n, const alpha_type alpha, const A_type* A, const OrdinalType lda, const B_type* B, const OrdinalType ldb, const beta_type beta, ScalarType* C, const OrdinalType ldc) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    alpha_type alpha_zero = ScalarTraits<alpha_type>::zero();
    beta_type beta_zero = ScalarTraits<beta_type>::zero();
    ScalarType C_zero = ScalarTraits<ScalarType>::zero();
    beta_type beta_one = ScalarTraits<beta_type>::one();
    OrdinalType i, j, k, NRowA = m;
    ScalarType temp1, temp2;
    bool BadArgument = false;
    bool Upper = (EUploChar[uplo] == 'U');
    if (ESideChar[side] == 'R') { NRowA = n; }

    // Quick return.
    if ( (m==izero) || (n==izero) || ( (alpha==alpha_zero)&&(beta==beta_one) ) ) { return; }
    if( m < izero ) {
      std::cout << "BLAS::SYMM Error: M == "<< m << std::endl;
      BadArgument = true; }
    if( n < izero ) {
      std::cout << "BLAS::SYMM Error: N == "<< n << std::endl;
      BadArgument = true; }
    if( lda < NRowA ) {
      std::cout << "BLAS::SYMM Error: LDA < "<<NRowA<<std::endl;
      BadArgument = true; }
    if( ldb < m ) {
      std::cout << "BLAS::SYMM Error: LDB < MAX(1,M)"<<std::endl;
      BadArgument = true; }
    if( ldc < m ) {
      std::cout << "BLAS::SYMM Error: LDC < MAX(1,M)"<<std::endl;
      BadArgument = true; }

    if(!BadArgument) {

      // Only need to scale C and return.
      if (alpha == alpha_zero) {
        if (beta == beta_zero ) {
          for (j=izero; j<n; j++) {
            for (i=izero; i<m; i++) {
              C[j*ldc + i] = C_zero;
            }
          }
        } else {
          for (j=izero; j<n; j++) {
            for (i=izero; i<m; i++) {
              C[j*ldc + i] *= beta;
            }
          }
        }
        return;
      }

      if ( ESideChar[side] == 'L') {
        // Compute C = alpha*A*B + beta*C

        if (Upper) {
          // The symmetric part of A is stored in the upper triangular part of the matrix.
          for (j=izero; j<n; j++) {
            for (i=izero; i<m; i++) {
              temp1 = alpha*B[j*ldb + i];
              temp2 = C_zero;
              for (k=izero; k<i; k++) {
                C[j*ldc + k] += temp1*A[i*lda + k];
                temp2 += B[j*ldb + k]*A[i*lda + k];
              }
              if (beta == beta_zero) {
                C[j*ldc + i] = temp1*A[i*lda + i] + alpha*temp2;
              } else {
                C[j*ldc + i] = beta*C[j*ldc + i] + temp1*A[i*lda + i] + alpha*temp2;
              }
            }
          }
        } else {
          // The symmetric part of A is stored in the lower triangular part of the matrix.
          for (j=izero; j<n; j++) {
            for (i=m-ione; i>-ione; i--) {
              temp1 = alpha*B[j*ldb + i];
              temp2 = C_zero;
              for (k=i+ione; k<m; k++) {
                C[j*ldc + k] += temp1*A[i*lda + k];
                temp2 += B[j*ldb + k]*A[i*lda + k];
              }
              if (beta == beta_zero) {
                C[j*ldc + i] = temp1*A[i*lda + i] + alpha*temp2;
              } else {
                C[j*ldc + i] = beta*C[j*ldc + i] + temp1*A[i*lda + i] + alpha*temp2;
              }
            }
          }
        }
      } else {
        // Compute C = alpha*B*A + beta*C.
        for (j=izero; j<n; j++) {
          temp1 = alpha*A[j*lda + j];
          if (beta == beta_zero) {
            for (i=izero; i<m; i++) {
              C[j*ldc + i] = temp1*B[j*ldb + i];
            }
          } else {
            for (i=izero; i<m; i++) {
              C[j*ldc + i] = beta*C[j*ldc + i] + temp1*B[j*ldb + i];
            }
          }
          for (k=izero; k<j; k++) {
            if (Upper) {
              temp1 = alpha*A[j*lda + k];
            } else {
              temp1 = alpha*A[k*lda + j];
            }
            for (i=izero; i<m; i++) {
              C[j*ldc + i] += temp1*B[k*ldb + i];
            }
          }
          for (k=j+ione; k<n; k++) {
            if (Upper) {
              temp1 = alpha*A[k*lda + j];
            } else {
              temp1 = alpha*A[j*lda + k];
            }
            for (i=izero; i<m; i++) {
              C[j*ldc + i] += temp1*B[k*ldb + i];
            }
          }
        }
      } // end if (ESideChar[side]=='L') ...
    } // end if(!BadArgument) ...
} // end SYMM

  template<typename OrdinalType, typename ScalarType>
  template <typename alpha_type, typename A_type, typename beta_type>
  void DefaultBLASImpl<OrdinalType, ScalarType>::SYRK(EUplo uplo, ETransp trans, const OrdinalType n, const OrdinalType k, const alpha_type alpha, const A_type* A, const OrdinalType lda, const beta_type beta, ScalarType* C, const OrdinalType ldc) const
  {
    typedef TypeNameTraits<OrdinalType> OTNT;
    typedef TypeNameTraits<ScalarType> STNT;

    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    alpha_type alpha_zero = ScalarTraits<alpha_type>::zero();
    beta_type beta_zero = ScalarTraits<beta_type>::zero();
    beta_type beta_one = ScalarTraits<beta_type>::one();
    A_type temp, A_zero = ScalarTraits<A_type>::zero();
    ScalarType C_zero = ScalarTraits<ScalarType>::zero();
    OrdinalType i, j, l, NRowA = n;
    bool BadArgument = false;
    bool Upper = (EUploChar[uplo] == 'U');

    TEUCHOS_TEST_FOR_EXCEPTION(
      Teuchos::ScalarTraits<ScalarType>::isComplex
      && (trans == CONJ_TRANS),
      std::logic_error,
            "Teuchos::BLAS<"<<OTNT::name()<<","<<STNT::name()<<">::SYRK()"
            " does not support CONJ_TRANS for complex data types."
      );

    // Change dimensions of A matrix is transposed
    if( !(ETranspChar[trans]=='N') ) {
      NRowA = k;
    }

    // Quick return.
    if ( n==izero ) { return; }
    if ( ( (alpha==alpha_zero) || (k==izero) )  && (beta==beta_one) ) { return; }
    if( n < izero ) {
      std::cout << "BLAS::SYRK Error: N == "<< n <<std::endl;
      BadArgument = true; }
    if( k < izero ) {
      std::cout << "BLAS::SYRK Error: K == "<< k <<std::endl;
      BadArgument = true; }
    if( lda < NRowA ) {
      std::cout << "BLAS::SYRK Error: LDA < "<<NRowA<<std::endl;
      BadArgument = true; }
    if( ldc < n ) {
      std::cout << "BLAS::SYRK Error: LDC < MAX(1,N)"<<std::endl;
      BadArgument = true; }

    if(!BadArgument) {

      // Scale C when alpha is zero
      if (alpha == alpha_zero) {
        if (Upper) {
          if (beta==beta_zero) {
            for (j=izero; j<n; j++) {
              for (i=izero; i<=j; i++) {
                C[j*ldc + i] = C_zero;
              }
            }
          }
          else {
            for (j=izero; j<n; j++) {
              for (i=izero; i<=j; i++) {
                C[j*ldc + i] *= beta;
              }
            }
          }
        }
        else {
          if (beta==beta_zero) {
            for (j=izero; j<n; j++) {
              for (i=j; i<n; i++) {
                C[j*ldc + i] = C_zero;
              }
            }
          }
          else {
            for (j=izero; j<n; j++) {
              for (i=j; i<n; i++) {
                C[j*ldc + i] *= beta;
              }
            }
          }
        }
        return;
      }

      // Now we can start the computation

      if ( ETranspChar[trans]=='N' ) {

        // Form C <- alpha*A*A' + beta*C
        if (Upper) {
          for (j=izero; j<n; j++) {
            if (beta==beta_zero) {
              for (i=izero; i<=j; i++) {
                C[j*ldc + i] = C_zero;
              }
            }
            else if (beta!=beta_one) {
              for (i=izero; i<=j; i++) {
                C[j*ldc + i] *= beta;
              }
            }
            for (l=izero; l<k; l++) {
              if (A[l*lda + j] != A_zero) {
                temp = alpha*A[l*lda + j];
                for (i = izero; i <=j; i++) {
                  C[j*ldc + i] += temp*A[l*lda + i];
                }
              }
            }
          }
        }
        else {
          for (j=izero; j<n; j++) {
            if (beta==beta_zero) {
              for (i=j; i<n; i++) {
                C[j*ldc + i] = C_zero;
              }
            }
            else if (beta!=beta_one) {
              for (i=j; i<n; i++) {
                C[j*ldc + i] *= beta;
              }
            }
            for (l=izero; l<k; l++) {
              if (A[l*lda + j] != A_zero) {
                temp = alpha*A[l*lda + j];
                for (i=j; i<n; i++) {
                  C[j*ldc + i] += temp*A[l*lda + i];
                }
              }
            }
          }
        }
      }
      else {

        // Form C <- alpha*A'*A + beta*C
        if (Upper) {
          for (j=izero; j<n; j++) {
            for (i=izero; i<=j; i++) {
              temp = A_zero;
              for (l=izero; l<k; l++) {
                temp += A[i*lda + l]*A[j*lda + l];
              }
              if (beta==beta_zero) {
                C[j*ldc + i] = alpha*temp;
              }
              else {
                C[j*ldc + i] = alpha*temp + beta*C[j*ldc + i];
              }
            }
          }
        }
        else {
          for (j=izero; j<n; j++) {
            for (i=j; i<n; i++) {
              temp = A_zero;
              for (l=izero; l<k; ++l) {
                temp += A[i*lda + l]*A[j*lda + l];
              }
              if (beta==beta_zero) {
                C[j*ldc + i] = alpha*temp;
              }
              else {
                C[j*ldc + i] = alpha*temp + beta*C[j*ldc + i];
              }
            }
          }
        }
      }
    } /* if (!BadArgument) */
  } /* END SYRK */

  template<typename OrdinalType, typename ScalarType>
  template <typename alpha_type, typename A_type>
  void DefaultBLASImpl<OrdinalType, ScalarType>::TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const alpha_type alpha, const A_type* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    alpha_type alpha_zero = ScalarTraits<alpha_type>::zero();
    A_type A_zero = ScalarTraits<A_type>::zero();
    ScalarType B_zero = ScalarTraits<ScalarType>::zero();
    ScalarType one = ScalarTraits<ScalarType>::one();
    OrdinalType i, j, k, NRowA = m;
    ScalarType temp;
    bool BadArgument = false;
    bool LSide = (ESideChar[side] == 'L');
    bool noUnit = (EDiagChar[diag] == 'N');
    bool Upper = (EUploChar[uplo] == 'U');
    bool noConj = (ETranspChar[transa] == 'T');

    if(!LSide) { NRowA = n; }

    // Quick return.
    if (n==izero || m==izero) { return; }
    if( m < izero ) {
      std::cout << "BLAS::TRMM Error: M == "<< m <<std::endl;
      BadArgument = true; }
    if( n < izero ) {
      std::cout << "BLAS::TRMM Error: N == "<< n <<std::endl;
      BadArgument = true; }
    if( lda < NRowA ) {
      std::cout << "BLAS::TRMM Error: LDA < "<<NRowA<<std::endl;
      BadArgument = true; }
    if( ldb < m ) {
      std::cout << "BLAS::TRMM Error: LDB < MAX(1,M)"<<std::endl;
      BadArgument = true; }

    if(!BadArgument) {

      // B only needs to be zeroed out.
      if( alpha == alpha_zero ) {
        for( j=izero; j<n; j++ ) {
          for( i=izero; i<m; i++ ) {
            B[j*ldb + i] = B_zero;
          }
        }
        return;
      }

      //  Start the computations.
      if ( LSide ) {
        // A is on the left side of B.

        if ( ETranspChar[transa]=='N' ) {
          // Compute B = alpha*A*B

          if ( Upper ) {
            // A is upper triangular
            for( j=izero; j<n; j++ ) {
              for( k=izero; k<m; k++) {
                if ( B[j*ldb + k] != B_zero ) {
                  temp = alpha*B[j*ldb + k];
                  for( i=izero; i<k; i++ ) {
                    B[j*ldb + i] += temp*A[k*lda + i];
                  }
                  if ( noUnit )
                    temp *=A[k*lda + k];
                  B[j*ldb + k] = temp;
                }
              }
            }
          } else {
            // A is lower triangular
            for( j=izero; j<n; j++ ) {
              for( k=m-ione; k>-ione; k-- ) {
                if( B[j*ldb + k] != B_zero ) {
                  temp = alpha*B[j*ldb + k];
                  B[j*ldb + k] = temp;
                  if ( noUnit )
                    B[j*ldb + k] *= A[k*lda + k];
                  for( i=k+ione; i<m; i++ ) {
                    B[j*ldb + i] += temp*A[k*lda + i];
                  }
                }
              }
            }
          }
        } else {
          // Compute B = alpha*A'*B or B = alpha*conj(A')*B
          if( Upper ) {
            for( j=izero; j<n; j++ ) {
              for( i=m-ione; i>-ione; i-- ) {
                temp = B[j*ldb + i];
                if ( noConj ) {
                  if( noUnit )
                    temp *= A[i*lda + i];
                  for( k=izero; k<i; k++ ) {
                    temp += A[i*lda + k]*B[j*ldb + k];
                  }
                } else {
                  if( noUnit )
                    temp *= ScalarTraits<A_type>::conjugate(A[i*lda + i]);
                  for( k=izero; k<i; k++ ) {
                    temp += ScalarTraits<A_type>::conjugate(A[i*lda + k])*B[j*ldb + k];
                  }
                }
                B[j*ldb + i] = alpha*temp;
              }
            }
          } else {
            for( j=izero; j<n; j++ ) {
              for( i=izero; i<m; i++ ) {
                temp = B[j*ldb + i];
                if ( noConj ) {
                  if( noUnit )
                    temp *= A[i*lda + i];
                  for( k=i+ione; k<m; k++ ) {
                    temp += A[i*lda + k]*B[j*ldb + k];
                  }
                } else {
                  if( noUnit )
                    temp *= ScalarTraits<A_type>::conjugate(A[i*lda + i]);
                  for( k=i+ione; k<m; k++ ) {
                    temp += ScalarTraits<A_type>::conjugate(A[i*lda + k])*B[j*ldb + k];
                  }
                }
                B[j*ldb + i] = alpha*temp;
              }
            }
          }
        }
      } else {
        // A is on the right hand side of B.

        if( ETranspChar[transa] == 'N' ) {
          // Compute B = alpha*B*A

          if( Upper ) {
            // A is upper triangular.
            for( j=n-ione; j>-ione; j-- ) {
              temp = alpha;
              if( noUnit )
                temp *= A[j*lda + j];
              for( i=izero; i<m; i++ ) {
                B[j*ldb + i] *= temp;
              }
              for( k=izero; k<j; k++ ) {
                if( A[j*lda + k] != A_zero ) {
                  temp = alpha*A[j*lda + k];
                  for( i=izero; i<m; i++ ) {
                    B[j*ldb + i] += temp*B[k*ldb + i];
                  }
                }
              }
            }
          } else {
            // A is lower triangular.
            for( j=izero; j<n; j++ ) {
              temp = alpha;
              if( noUnit )
                temp *= A[j*lda + j];
              for( i=izero; i<m; i++ ) {
                B[j*ldb + i] *= temp;
              }
              for( k=j+ione; k<n; k++ ) {
                if( A[j*lda + k] != A_zero ) {
                  temp = alpha*A[j*lda + k];
                  for( i=izero; i<m; i++ ) {
                    B[j*ldb + i] += temp*B[k*ldb + i];
                  }
                }
              }
            }
          }
        } else {
          // Compute B = alpha*B*A' or B = alpha*B*conj(A')

          if( Upper ) {
            for( k=izero; k<n; k++ ) {
              for( j=izero; j<k; j++ ) {
                if( A[k*lda + j] != A_zero ) {
                  if ( noConj )
                    temp = alpha*A[k*lda + j];
                  else
                    temp = alpha*ScalarTraits<A_type>::conjugate(A[k*lda + j]);
                  for( i=izero; i<m; i++ ) {
                    B[j*ldb + i] += temp*B[k*ldb + i];
                  }
                }
              }
              temp = alpha;
              if( noUnit ) {
                if ( noConj )
                  temp *= A[k*lda + k];
                else
                  temp *= ScalarTraits<A_type>::conjugate(A[k*lda + k]);
              }
              if( temp != one ) {
                for( i=izero; i<m; i++ ) {
                  B[k*ldb + i] *= temp;
                }
              }
            }
          } else {
            for( k=n-ione; k>-ione; k-- ) {
              for( j=k+ione; j<n; j++ ) {
                if( A[k*lda + j] != A_zero ) {
                  if ( noConj )
                    temp = alpha*A[k*lda + j];
                  else
                    temp = alpha*ScalarTraits<A_type>::conjugate(A[k*lda + j]);
                  for( i=izero; i<m; i++ ) {
                    B[j*ldb + i] += temp*B[k*ldb + i];
                  }
                }
              }
              temp = alpha;
              if( noUnit ) {
                if ( noConj )
                  temp *= A[k*lda + k];
                else
                  temp *= ScalarTraits<A_type>::conjugate(A[k*lda + k]);
              }
              if( temp != one ) {
                for( i=izero; i<m; i++ ) {
                  B[k*ldb + i] *= temp;
                }
              }
            }
          }
        } // end if( ETranspChar[transa] == 'N' ) ...
      } // end if ( LSide ) ...
    } // end if (!BadArgument)
  } // end TRMM

  template<typename OrdinalType, typename ScalarType>
  template <typename alpha_type, typename A_type>
  void DefaultBLASImpl<OrdinalType, ScalarType>::TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType m, const OrdinalType n, const alpha_type alpha, const A_type* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    alpha_type alpha_zero = ScalarTraits<alpha_type>::zero();
    A_type A_zero = ScalarTraits<A_type>::zero();
    ScalarType B_zero = ScalarTraits<ScalarType>::zero();
    alpha_type alpha_one = ScalarTraits<alpha_type>::one();
    ScalarType B_one = ScalarTraits<ScalarType>::one();
    ScalarType temp;
    OrdinalType NRowA = m;
    bool BadArgument = false;
    bool noUnit = (EDiagChar[diag]=='N');
    bool noConj = (ETranspChar[transa] == 'T');

    if (!(ESideChar[side] == 'L')) { NRowA = n; }

    // Quick return.
    if (n == izero || m == izero) { return; }
    if( m < izero ) {
      std::cout << "BLAS::TRSM Error: M == "<<m<<std::endl;
      BadArgument = true; }
    if( n < izero ) {
      std::cout << "BLAS::TRSM Error: N == "<<n<<std::endl;
      BadArgument = true; }
    if( lda < NRowA ) {
      std::cout << "BLAS::TRSM Error: LDA < "<<NRowA<<std::endl;
      BadArgument = true; }
    if( ldb < m ) {
      std::cout << "BLAS::TRSM Error: LDB < MAX(1,M)"<<std::endl;
      BadArgument = true; }

    if(!BadArgument)
      {
        int i, j, k;
        // Set the solution to the zero vector.
        if(alpha == alpha_zero) {
            for(j = izero; j < n; j++) {
                for( i = izero; i < m; i++) {
                    B[j*ldb+i] = B_zero;
                }
            }
        }
        else
        { // Start the operations.
            if(ESideChar[side] == 'L') {
                //
                // Perform computations for OP(A)*X = alpha*B
                //
                if(ETranspChar[transa] == 'N') {
                    //
                    //  Compute B = alpha*inv( A )*B
                    //
                    if(EUploChar[uplo] == 'U') {
                        // A is upper triangular.
                        for(j = izero; j < n; j++) {
                            // Perform alpha*B if alpha is not 1.
                          if(alpha != alpha_one) {
                                for( i = izero; i < m; i++) {
                                    B[j*ldb+i] *= alpha;
                                }
                            }
                            // Perform a backsolve for column j of B.
                            for(k = (m - ione); k > -ione; k--) {
                                // If this entry is zero, we don't have to do anything.
                                if (B[j*ldb + k] != B_zero) {
                                    if ( noUnit ) {
                                        B[j*ldb + k] /= A[k*lda + k];
                                    }
                                    for(i = izero; i < k; i++) {
                                        B[j*ldb + i] -= B[j*ldb + k] * A[k*lda + i];
                                    }
                                }
                            }
                        }
                    }
                    else
                    { // A is lower triangular.
                        for(j = izero; j < n; j++) {
                            // Perform alpha*B if alpha is not 1.
                            if(alpha != alpha_one) {
                                for( i = izero; i < m; i++) {
                                    B[j*ldb+i] *= alpha;
                                }
                            }
                            // Perform a forward solve for column j of B.
                            for(k = izero; k < m; k++) {
                                // If this entry is zero, we don't have to do anything.
                                if (B[j*ldb + k] != B_zero) {
                                    if ( noUnit ) {
                                        B[j*ldb + k] /= A[k*lda + k];
                                    }
                                    for(i = k+ione; i < m; i++) {
                                        B[j*ldb + i] -= B[j*ldb + k] * A[k*lda + i];
                                    }
                                }
                            }
                        }
                    } // end if (uplo == 'U')
                }  // if (transa =='N')
                else {
                    //
                    //  Compute B = alpha*inv( A' )*B
                    //  or      B = alpha*inv( conj(A') )*B
                    //
                    if(EUploChar[uplo] == 'U') {
                        // A is upper triangular.
                        for(j = izero; j < n; j++) {
                            for( i = izero; i < m; i++) {
                                temp = alpha*B[j*ldb+i];
                                if ( noConj ) {
                                  for(k = izero; k < i; k++) {
                                      temp -= A[i*lda + k] * B[j*ldb + k];
                                  }
                                  if ( noUnit ) {
                                      temp /= A[i*lda + i];
                                  }
                                } else {
                                  for(k = izero; k < i; k++) {
                                      temp -= ScalarTraits<A_type>::conjugate(A[i*lda + k])
                                            * B[j*ldb + k];
                                  }
                                  if ( noUnit ) {
                                      temp /= ScalarTraits<A_type>::conjugate(A[i*lda + i]);
                                  }
                                }
                                B[j*ldb + i] = temp;
                            }
                        }
                    }
                    else
                    { // A is lower triangular.
                        for(j = izero; j < n; j++) {
                            for(i = (m - ione) ; i > -ione; i--) {
                                temp = alpha*B[j*ldb+i];
                                if ( noConj ) {
                                  for(k = i+ione; k < m; k++) {
                                    temp -= A[i*lda + k] * B[j*ldb + k];
                                  }
                                  if ( noUnit ) {
                                    temp /= A[i*lda + i];
                                  }
                                } else {
                                  for(k = i+ione; k < m; k++) {
                                    temp -= ScalarTraits<A_type>::conjugate(A[i*lda + k])
                                          * B[j*ldb + k];
                                  }
                                  if ( noUnit ) {
                                    temp /= ScalarTraits<A_type>::conjugate(A[i*lda + i]);
                                  }
                                }
                                B[j*ldb + i] = temp;
                            }
                        }
                    }
                }
            }  // if (side == 'L')
            else {
               // side == 'R'
               //
               // Perform computations for X*OP(A) = alpha*B
               //
              if (ETranspChar[transa] == 'N') {
                    //
                    //  Compute B = alpha*B*inv( A )
                    //
                    if(EUploChar[uplo] == 'U') {
                        // A is upper triangular.
                        // Perform a backsolve for column j of B.
                        for(j = izero; j < n; j++) {
                            // Perform alpha*B if alpha is not 1.
                            if(alpha != alpha_one) {
                                for( i = izero; i < m; i++) {
                                    B[j*ldb+i] *= alpha;
                                }
                            }
                            for(k = izero; k < j; k++) {
                                // If this entry is zero, we don't have to do anything.
                                if (A[j*lda + k] != A_zero) {
                                    for(i = izero; i < m; i++) {
                                        B[j*ldb + i] -= A[j*lda + k] * B[k*ldb + i];
                                    }
                                }
                            }
                            if ( noUnit ) {
                                temp = B_one/A[j*lda + j];
                                for(i = izero; i < m; i++) {
                                    B[j*ldb + i] *= temp;
                                }
                            }
                        }
                    }
                    else
                    { // A is lower triangular.
                        for(j = (n - ione); j > -ione; j--) {
                            // Perform alpha*B if alpha is not 1.
                            if(alpha != alpha_one) {
                                for( i = izero; i < m; i++) {
                                    B[j*ldb+i] *= alpha;
                                }
                            }
                            // Perform a forward solve for column j of B.
                            for(k = j+ione; k < n; k++) {
                                // If this entry is zero, we don't have to do anything.
                                if (A[j*lda + k] != A_zero) {
                                    for(i = izero; i < m; i++) {
                                        B[j*ldb + i] -= A[j*lda + k] * B[k*ldb + i];
                                    }
                                }
                            }
                            if ( noUnit ) {
                                temp = B_one/A[j*lda + j];
                                for(i = izero; i < m; i++) {
                                    B[j*ldb + i] *= temp;
                                }
                            }
                        }
                    } // end if (uplo == 'U')
                }  // if (transa =='N')
                else {
                    //
                    //  Compute B = alpha*B*inv( A' )
                    //  or      B = alpha*B*inv( conj(A') )
                    //
                    if(EUploChar[uplo] == 'U') {
                        // A is upper triangular.
                        for(k = (n - ione); k > -ione; k--) {
                            if ( noUnit ) {
                                if ( noConj )
                                  temp = B_one/A[k*lda + k];
                                else
                                  temp = B_one/ScalarTraits<A_type>::conjugate(A[k*lda + k]);
                                for(i = izero; i < m; i++) {
                                    B[k*ldb + i] *= temp;
                                }
                            }
                            for(j = izero; j < k; j++) {
                                if (A[k*lda + j] != A_zero) {
                                    if ( noConj )
                                      temp = A[k*lda + j];
                                    else
                                      temp = ScalarTraits<A_type>::conjugate(A[k*lda + j]);
                                    for(i = izero; i < m; i++) {
                                        B[j*ldb + i] -= temp*B[k*ldb + i];
                                    }
                                }
                            }
                            if (alpha != alpha_one) {
                                for (i = izero; i < m; i++) {
                                    B[k*ldb + i] *= alpha;
                                }
                            }
                        }
                    }
                    else
                    { // A is lower triangular.
                        for(k = izero; k < n; k++) {
                            if ( noUnit ) {
                                if ( noConj )
                                  temp = B_one/A[k*lda + k];
                                else
                                  temp = B_one/ScalarTraits<A_type>::conjugate(A[k*lda + k]);
                                for (i = izero; i < m; i++) {
                                    B[k*ldb + i] *= temp;
                                }
                            }
                            for(j = k+ione; j < n; j++) {
                                if(A[k*lda + j] != A_zero) {
                                    if ( noConj )
                                      temp = A[k*lda + j];
                                    else
                                      temp = ScalarTraits<A_type>::conjugate(A[k*lda + j]);
                                    for(i = izero; i < m; i++) {
                                        B[j*ldb + i] -= temp*B[k*ldb + i];
                                    }
                                }
                            }
                            if (alpha != alpha_one) {
                                for (i = izero; i < m; i++) {
                                    B[k*ldb + i] *= alpha;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
  }

  // Explicit instantiation for template<int,float>

  template <>
  class TEUCHOSNUMERICS_LIB_DLL_EXPORT BLAS<int, float>
  {
  public:
    inline BLAS(void) {}
    inline BLAS(const BLAS<int, float>& /*BLAS_source*/) {}
    inline virtual ~BLAS(void) {}
    void ROTG(float* da, float* db, float* c, float* s) const;
    void ROT(const int n, float* dx, const int incx, float* dy, const int incy, float* c, float* s) const;
    float ASUM(const int n, const float* x, const int incx) const;
    void AXPY(const int n, const float alpha, const float* x, const int incx, float* y, const int incy) const;
    void COPY(const int n, const float* x, const int incx, float* y, const int incy) const;
    float DOT(const int n, const float* x, const int incx, const float* y, const int incy) const;
    float NRM2(const int n, const float* x, const int incx) const;
    void SCAL(const int n, const float alpha, float* x, const int incx) const;
    int IAMAX(const int n, const float* x, const int incx) const;
    void GEMV(ETransp trans, const int m, const int n, const float alpha, const float* A, const int lda, const float* x, const int incx, const float beta, float* y, const int incy) const;
    void TRMV(EUplo uplo, ETransp trans, EDiag diag, const int n, const float* A, const int lda, float* x, const int incx) const;
    void GER(const int m, const int n, const float alpha, const float* x, const int incx, const float* y, const int incy, float* A, const int lda) const;
    void GEMM(ETransp transa, ETransp transb, const int m, const int n, const int k, const float alpha, const float* A, const int lda, const float* B, const int ldb, const float beta, float* C, const int ldc) const;
    void SYMM(ESide side, EUplo uplo, const int m, const int n, const float alpha, const float* A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc) const;
    void SYRK(EUplo uplo, ETransp trans, const int n, const int k, const float alpha, const float* A, const int lda, const float beta, float* C, const int ldc) const;
    void TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const float alpha, const float* A, const int lda, float* B, const int ldb) const;
    void TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const float alpha, const float* A, const int lda, float* B, const int ldb) const;
  };

  // Explicit instantiation for template<int,double>

  template<>
  class TEUCHOSNUMERICS_LIB_DLL_EXPORT BLAS<int, double>
  {
  public:
    inline BLAS(void) {}
    inline BLAS(const BLAS<int, double>& /*BLAS_source*/) {}
    inline virtual ~BLAS(void) {}
    void ROTG(double* da, double* db, double* c, double* s) const;
    void ROT(const int n, double* dx, const int incx, double* dy, const int incy, double* c, double* s) const;
    double ASUM(const int n, const double* x, const int incx) const;
    void AXPY(const int n, const double alpha, const double* x, const int incx, double* y, const int incy) const;
    void COPY(const int n, const double* x, const int incx, double* y, const int incy) const;
    double DOT(const int n, const double* x, const int incx, const double* y, const int incy) const;
    double NRM2(const int n, const double* x, const int incx) const;
    void SCAL(const int n, const double alpha, double* x, const int incx) const;
    int IAMAX(const int n, const double* x, const int incx) const;
    void GEMV(ETransp trans, const int m, const int n, const double alpha, const double* A, const int lda, const double* x, const int incx, const double beta, double* y, const int incy) const;
    void TRMV(EUplo uplo, ETransp trans, EDiag diag, const int n, const double* A, const int lda, double* x, const int incx) const;
    void GER(const int m, const int n, const double alpha, const double* x, const int incx, const double* y, const int incy, double* A, const int lda) const;
    void GEMM(ETransp transa, ETransp transb, const int m, const int n, const int k, const double alpha, const double* A, const int lda, const double* B, const int ldb, const double beta, double* C, const int ldc) const;
    void SYMM(ESide side, EUplo uplo, const int m, const int n, const double alpha, const double* A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) const;
    void SYRK(EUplo uplo, ETransp trans, const int n, const int k, const double alpha, const double* A, const int lda, const double beta, double* C, const int ldc) const;
    void TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const double alpha, const double* A, const int lda, double* B, const int ldb) const;
    void TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const double alpha, const double* A, const int lda, double* B, const int ldb) const;
  };

  // Explicit instantiation for template<int,complex<float> >

  template<>
  class TEUCHOSNUMERICS_LIB_DLL_EXPORT BLAS<int, std::complex<float> >
  {
  public:
    inline BLAS(void) {}
    inline BLAS(const BLAS<int, std::complex<float> >& /*BLAS_source*/) {}
    inline virtual ~BLAS(void) {}
    void ROTG(std::complex<float>* da, std::complex<float>* db, float* c, std::complex<float>* s) const;
    void ROT(const int n, std::complex<float>* dx, const int incx, std::complex<float>* dy, const int incy, float* c, std::complex<float>* s) const;
    float ASUM(const int n, const std::complex<float>* x, const int incx) const;
    void AXPY(const int n, const std::complex<float> alpha, const std::complex<float>* x, const int incx, std::complex<float>* y, const int incy) const;
    void COPY(const int n, const std::complex<float>* x, const int incx, std::complex<float>* y, const int incy) const;
    std::complex<float> DOT(const int n, const std::complex<float>* x, const int incx, const std::complex<float>* y, const int incy) const;
    float NRM2(const int n, const std::complex<float>* x, const int incx) const;
    void SCAL(const int n, const std::complex<float> alpha, std::complex<float>* x, const int incx) const;
    int IAMAX(const int n, const std::complex<float>* x, const int incx) const;
    void GEMV(ETransp trans, const int m, const int n, const std::complex<float> alpha, const std::complex<float>* A, const int lda, const std::complex<float>* x, const int incx, const std::complex<float> beta, std::complex<float>* y, const int incy) const;
    void TRMV(EUplo uplo, ETransp trans, EDiag diag, const int n, const std::complex<float>* A, const int lda, std::complex<float>* x, const int incx) const;
    void GER(const int m, const int n, const std::complex<float> alpha, const std::complex<float>* x, const int incx, const std::complex<float>* y, const int incy, std::complex<float>* A, const int lda) const;
    void GEMM(ETransp transa, ETransp transb, const int m, const int n, const int k, const std::complex<float> alpha, const std::complex<float>* A, const int lda, const std::complex<float>* B, const int ldb, const std::complex<float> beta, std::complex<float>* C, const int ldc) const;
    void SYMM(ESide side, EUplo uplo, const int m, const int n, const std::complex<float> alpha, const std::complex<float>* A, const int lda, const std::complex<float> *B, const int ldb, const std::complex<float> beta, std::complex<float> *C, const int ldc) const;
    void SYRK(EUplo uplo, ETransp trans, const int n, const int k, const std::complex<float> alpha, const std::complex<float>* A, const int lda, const std::complex<float> beta, std::complex<float>* C, const int ldc) const;
    void TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const std::complex<float> alpha, const std::complex<float>* A, const int lda, std::complex<float>* B, const int ldb) const;
    void TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const std::complex<float> alpha, const std::complex<float>* A, const int lda, std::complex<float>* B, const int ldb) const;
  };

  // Explicit instantiation for template<int,complex<double> >

  template<>
  class TEUCHOSNUMERICS_LIB_DLL_EXPORT BLAS<int, std::complex<double> >
  {
  public:
    inline BLAS(void) {}
    inline BLAS(const BLAS<int, std::complex<double> >& /*BLAS_source*/) {}
    inline virtual ~BLAS(void) {}
    void ROTG(std::complex<double>* da, std::complex<double>* db, double* c, std::complex<double>* s) const;
    void ROT(const int n, std::complex<double>* dx, const int incx, std::complex<double>* dy, const int incy, double* c, std::complex<double>* s) const;
    double ASUM(const int n, const std::complex<double>* x, const int incx) const;
    void AXPY(const int n, const std::complex<double> alpha, const std::complex<double>* x, const int incx, std::complex<double>* y, const int incy) const;
    void COPY(const int n, const std::complex<double>* x, const int incx, std::complex<double>* y, const int incy) const;
    std::complex<double> DOT(const int n, const std::complex<double>* x, const int incx, const std::complex<double>* y, const int incy) const;
    double NRM2(const int n, const std::complex<double>* x, const int incx) const;
    void SCAL(const int n, const std::complex<double> alpha, std::complex<double>* x, const int incx) const;
    int IAMAX(const int n, const std::complex<double>* x, const int incx) const;
    void GEMV(ETransp trans, const int m, const int n, const std::complex<double> alpha, const std::complex<double>* A, const int lda, const std::complex<double>* x, const int incx, const std::complex<double> beta, std::complex<double>* y, const int incy) const;
    void TRMV(EUplo uplo, ETransp trans, EDiag diag, const int n, const std::complex<double>* A, const int lda, std::complex<double>* x, const int incx) const;
    void GER(const int m, const int n, const std::complex<double> alpha, const std::complex<double>* x, const int incx, const std::complex<double>* y, const int incy, std::complex<double>* A, const int lda) const;
    void GEMM(ETransp transa, ETransp transb, const int m, const int n, const int k, const std::complex<double> alpha, const std::complex<double>* A, const int lda, const std::complex<double>* B, const int ldb, const std::complex<double> beta, std::complex<double>* C, const int ldc) const;
    void SYMM(ESide side, EUplo uplo, const int m, const int n, const std::complex<double> alpha, const std::complex<double>* A, const int lda, const std::complex<double> *B, const int ldb, const std::complex<double> beta, std::complex<double> *C, const int ldc) const;
    void SYRK(EUplo uplo, ETransp trans, const int n, const int k, const std::complex<double> alpha, const std::complex<double>* A, const int lda, const std::complex<double> beta, std::complex<double>* C, const int ldc) const;
    void TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const std::complex<double> alpha, const std::complex<double>* A, const int lda, std::complex<double>* B, const int ldb) const;
    void TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int m, const int n, const std::complex<double> alpha, const std::complex<double>* A, const int lda, std::complex<double>* B, const int ldb) const;
  };

} // namespace Teuchos

#endif // _TEUCHOS_BLAS_HPP_
