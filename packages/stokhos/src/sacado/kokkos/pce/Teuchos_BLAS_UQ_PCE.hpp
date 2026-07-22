// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _TEUCHOS_BLAS_UQ_PCE_HPP_
#define _TEUCHOS_BLAS_UQ_PCE_HPP_

#include "Teuchos_BLAS.hpp"
#include "Sacado_UQ_PCE.hpp"

// Specialize some things used in the default BLAS implementation that
// don't seem correct for UQ::PCE scalar type
namespace Teuchos {

  namespace details {

    template<typename Storage>
    class GivensRotator<Sacado::UQ::PCE<Storage>, false> {
    public:
      typedef Sacado::UQ::PCE<Storage> ScalarType;
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

#endif // _TEUCHOS_BLAS_UQ_PCE_HPP_
