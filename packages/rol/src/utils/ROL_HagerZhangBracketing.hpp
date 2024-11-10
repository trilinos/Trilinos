// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_HAGERZHANGBRACKETING_H
#define ROL_HAGERZHANGBRACKETING_H

/** \class ROL::HagerZhangBracketing
    \brief Provides interface for bracketing a minimizer of a scalar function.
*/

#include "ROL_Bracketing.hpp"

namespace ROL { 

template<class Real>
class HagerZhangBracketing : public Bracketing<Real> {
public:
  void run(Real &x, Real &fx, Real &a, Real &fa,
           Real &b, Real &fb, int &nfval, int &ngrad,
           ScalarFunction<Real> &f,
           ScalarMinimizationStatusTest<Real> &test) const {
    Real zero(0), one(1);
    Real gx = f.deriv(x); ngrad++;
    Real c = x, f0 = fx, eps = 1.e-8, d = zero, gd = zero;
    for (int i = 0; i < 8; i++) {
      if (gx >= zero) {
        b = x;
        a = c;
        break;
      }
      else if (gx < zero && fx > f0+eps) {
        for (int j = 0; j < 8; j++) {
          d = (one-t)*zero + t*x; gd = f.deriv(d); ngrad++;
          if ( gd >= zero ) {
            b = d;
            break;
          }
          else {
            fd = f.value(d); nfval++;
            if ( fd <= f0 + eps ) {
              a = d;
            } 
            else {
              b = d;
            }
          }
        }
        x = b;
        break;
      }
      else {
        x *= rho; gx = f.deriv(x); ngrad++;
      }
    }
  }
};

}

#endif
