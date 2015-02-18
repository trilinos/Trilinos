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

#ifndef _TEUCHOS_BLAS_MP_VECTOR_HPP_
#define _TEUCHOS_BLAS_MP_VECTOR_HPP_

#include "Teuchos_BLAS.hpp"
#include "Sacado_MP_Vector.hpp"

// Specialize some things used in the default BLAS implementation that
// don't seem correct for MP::Vector scalar type
namespace Teuchos {

  namespace details {

    template<typename Storage>
    class GivensRotator<Sacado::MP::Vector<Storage>, false> {
    public:
      typedef Sacado::MP::Vector<Storage> ScalarType;
      typedef ScalarType c_type;

      void
      ROTG (ScalarType* da,
            ScalarType* db,
            ScalarType* c,
            ScalarType* s) const {

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

    private:

      /// Return ABS(x) if y > 0 or y is +0, else -ABS(x) (if y is -0 or < 0).
      ScalarType SIGN (ScalarType x, ScalarType y) const {
        typedef typename ScalarType::value_type value_type;
        typedef typename ScalarType::ordinal_type ordinal_type;

        GivensRotator<value_type> value_rotator;
        const ordinal_type sz = x.size() > y.size() ? x.size() : y.size();
        ScalarType z(sz, 0.0);
        for (ordinal_type i=0; i<sz; ++i)
          z.fastAccessCoeff(i) = value_rotator.SIGN(x.coeff(i), y.coeff(i));
        return z;
      }
    };

  } // namespace details

} // namespace Teuchos

#endif // _TEUCHOS_BLAS__MP_VECTOR_HPP_
