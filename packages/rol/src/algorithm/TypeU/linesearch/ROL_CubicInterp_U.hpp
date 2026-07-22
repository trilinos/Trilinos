// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CUBICINTERP_U_H
#define ROL_CUBICINTERP_U_H

/** \class ROL::CubicInterp_U
    \brief Implements cubic interpolation back tracking line search.
*/

#include "ROL_LineSearch_U.hpp"

namespace ROL {

template<typename Real>
class CubicInterp_U : public LineSearch_U<Real> {
private:
  Real rho_;
  Ptr<Vector<Real>> xnew_; 

  using LineSearch_U<Real>::getInitialAlpha;
  using LineSearch_U<Real>::status;

public:

  CubicInterp_U( ParameterList &parlist ) : LineSearch_U<Real>(parlist) { 
    const Real half(0.5);
    rho_ = parlist.sublist("Step").sublist("Line Search").sublist("Line-Search Method").get("Backtracking Rate",half);
  }

  void initialize(const Vector<Real> &x, const Vector<Real> &g) override {
    LineSearch_U<Real>::initialize(x,g);
    xnew_   = x.clone();
  }

  void run( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
            const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
            Objective<Real> &obj ) override {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    ls_neval = 0;
    ls_ngrad = 0;
    // Get initial line search parameter
    alpha = getInitialAlpha(ls_neval,ls_ngrad,fval,gs,x,s,obj);
    // Update iterate
    xnew_->set(x); xnew_->axpy(alpha,s);
    // Get objective value at xnew
    Real fold = fval;
    obj.update(*xnew_,UpdateType::Trial);
    fval = obj.value(*xnew_,tol);
    ls_neval++;
    // Initialize
    Real fvalp(0), alpha1(0), alpha2(0), a(0), b(0), x1(0), x2(0);
    const Real one(1), two(2), three(3), half(0.5), p1(0.1);
    bool first_iter = true;
    // Perform cubic interpolation back tracking
    while (!status(LINESEARCH_U_CUBICINTERP,ls_neval,ls_ngrad,alpha,fold,gs,fval,x,s,obj)) {
      if (first_iter) { // Minimize quadratic interpolate
        alpha1 = -gs*alpha*alpha/(two*(fval-fold-gs*alpha));
        first_iter = false;
      }
      else {              // Minimize cubic interpolate
        x1 = fval-fold-alpha*gs;
        x2 = fvalp-fval-alpha2*gs;
        a = (one/(alpha - alpha2))*( x1/(alpha*alpha) - x2/(alpha2*alpha2));
        b = (one/(alpha - alpha2))*(-x1*alpha2/(alpha*alpha) + x2*alpha/(alpha2*alpha2));
        if (std::abs(a) < ROL_EPSILON<Real>()) {
          alpha1 = -gs/(two*b);
        }
        else {
          alpha1 = (-b+sqrt(b*b-three*a*gs))/(three*a);
        }
        if ( alpha1 > half*alpha ) {
          alpha1 = half*alpha;
        }
      }
      alpha2 = alpha;
      fvalp  = fval;
      // Back track if necessary
      if (alpha1 <= p1*alpha) {
        alpha *= p1;
      }
      else if (alpha1 >= half*alpha) {
        alpha *= half;
      }
      else {
        alpha = alpha1;
      }
      // Update iterate
      xnew_->set(x); xnew_->axpy(alpha,s);
      // Get objective value at xnew
      obj.update(*xnew_,UpdateType::Trial);
      fval = obj.value(*xnew_,tol);
      ls_neval++;
    }
  }
}; // class ROL::CubicInterp_U

} // namespace ROL

#endif
