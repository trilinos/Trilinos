// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_GOLDENSECTIONSCALARMINIMIZATION_H
#define ROL_GOLDENSECTIONSCALARMINIMIZATION_H

/** \class ROL::GoldenSectionScalarMinimization
    \brief Implements the golden section method to minimize a scalar function on a
           bounded interval.
*/

#include "ROL_ScalarMinimization.hpp"
#include "ROL_ScalarFunction.hpp"
#include "ROL_Types.hpp"

namespace ROL { 

template<class Real>
class GoldenSectionScalarMinimization : public ScalarMinimization<Real> {
private:
  Real tol_;
  int niter_;

public:
  // Constructor
  GoldenSectionScalarMinimization( ROL::ParameterList &parlist ) {
    ROL::ParameterList &list = parlist.sublist("Scalar Minimization").sublist("Golden Section");
    tol_   = list.get("Tolerance",1.e-10);
    niter_ = list.get("Iteration Limit",1000);
  }

  void run(Real &fx, Real &x, int &nfval, int &ngrad,
           ScalarFunction<Real> &f,
           const Real A, const Real B,
           ScalarMinimizationStatusTest<Real> &test) const {
    Real one(1), two(2), five(5);
    nfval = 0; ngrad = 0;
    // Reciprocal of golden ratio
    const Real c = two/(one + std::sqrt(five));
    // Compute value f(a), f(b), f(u), and f(v)
    Real a = A,               fa = f.value(a); nfval++;
    Real b = B,               fb = f.value(b); nfval++;
    Real u = c*a + (one-c)*b, fu = f.value(u); nfval++;
    Real v = (one-c)*a + c*b, fv = f.value(v); nfval++;
    Real gx = ROL_INF<Real>();
    bool deriv = false;
    // Find minimum of all function evaluations
    if ( fa <= fu && fa <= fv && fa <= fb ) {
      x = a; fx = fa;
    }
    else if ( fu <= fa && fu <= fv && fu <= fb ) {
      x = u; fx = fu;
    }
    else if ( fv <= fa && fv <= fu && fv <= fb ) {
      x = v; fx = fv;
    }
    else {
      x = b; fx = fb;
    }
    // Run Golden Section
    for ( int i = 0; i < niter_; i++ ) {
      if ( std::abs(b - a) < tol_ || test.check(x,fx,gx,nfval,ngrad,deriv) ) {
        break;
      }
      if ( fu > fv ) {
        a = u;               fa = fu;
        u = v;               fu = fv;
        v = (one-c)*a + c*b; fv = f.value(v); nfval++;
      }
      else {
        b = v;               fb = fv;
        v = u;               fv = fu;
        u = c*a + (one-c)*b; fu = f.value(u); nfval++;
      }
      if ( fa <= fu && fa <= fv && fa <= fb ) {
        x = a; fx = fa;
      }
      else if ( fu <= fa && fu <= fv && fu <= fb ) {
        x = u; fx = fu;
      }
      else if ( fv <= fa && fv <= fu && fv <= fb ) {
        x = v; fx = fv;
      }
      else {
        x = b; fx = fb;
      }
    }
  }
};

}

#endif
