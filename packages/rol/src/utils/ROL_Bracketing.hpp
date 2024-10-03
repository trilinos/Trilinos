// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BRACKETING_H
#define ROL_BRACKETING_H

/** \class ROL::Bracketing
    \brief Provides interface for bracketing a minimizer of a scalar function.
*/

#include "ROL_ScalarFunction.hpp"
#include "ROL_ScalarMinimizationStatusTest.hpp"

namespace ROL { 

template<class Real>
class Bracketing {
public:
  virtual ~Bracketing() {}
  void run(Real &x, Real &fx, Real &a, Real &fa,
           Real &b, Real &fb, int &nfval, int &ngrad,
           ScalarFunction<Real> &f) const {
    ScalarMinimizationStatusTest<Real> test;
    run(x,fx,a,fa,b,fb,nfval,ngrad,f,test);
  }
  virtual void run(Real &x, Real &fx, Real &a, Real &fa,
                   Real &b, Real &fb, int &nfval, int &ngrad,
                   ScalarFunction<Real> &f,
                   ScalarMinimizationStatusTest<Real> &test) const {
    Real zero(0), half(0.5), one(1);
    const Real c(1.618034);
    const Real eps(ROL_EPSILON<Real>());
    const Real lim(100);
    Real r = zero, q = zero, u = zero, v = zero, ulim = zero, fu = zero, gx = ROL_INF<Real>();
    bool deriv = false;
    // f(a) is assumed to be greater than or equal to f(b)
    if ( fb > fa ) {
      return;
    }
    x = b + c*(b-a); fx = f.value(x); nfval++;
    for ( int i = 0; i < 8; i++ ) {
      if (fb < fx || test.check(x,fx,gx,nfval,ngrad,deriv) ) {
        break;
      }
      r = (b-a)*(fb-fx);
      q = (b-x)*(fb-fa);
      v = ((q > r) ? one : -one)*std::max(std::abs(q-r),eps);
      u = b - half*((b-x)*q - (b-a)*r)/v;
      ulim = b + lim*(x-b);
      if ( (b-u)*(u-x) > zero ) {
        fu = f.value(u); nfval++;
        if ( fu < fx ) {
          a = b; fa = fb;
          b = u; fb = fu;
          break;
        }
        else if ( fu > fb ) {
          x = u; fx = fu;
          break;
        }
        u = x + c*(x-b); fu = f.value(u); nfval++;
      }
      else if ( (x-u)*(u-ulim) > zero ) {
        fu = f.value(u); nfval++;
        if ( fu < fx ) {
          b = x; fb = fx;
          x = u; fx = fu;
          u = x + c*(x-b); fu = f.value(u); nfval++;
        }
      }
      else if ( (u-ulim)*(ulim-x) > zero) {
        u = ulim; fu = f.value(u); nfval++;
      }
      else {
        u = x + c*(x-b); fu = f.value(u); nfval++;
      }
      a = b; fa = fb;
      b = x; fb = fx;
      x = u; fx = fu;
    }
  }
};

}

#endif
