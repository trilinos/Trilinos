// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BACKTRACKING_U_H
#define ROL_BACKTRACKING_U_H

#include "ROL_LineSearch_U.hpp"

/** \class ROL::BackTracking_U
    \brief Implements a simple back tracking line search.
*/

namespace ROL {

template<typename Real>
class BackTracking_U : public LineSearch_U<Real> {
private:
  Real rho_;
  Ptr<Vector<Real>> xnew_; 

  using LineSearch_U<Real>::getInitialAlpha;
  using LineSearch_U<Real>::status;

public:

  BackTracking_U(ParameterList &parlist) : LineSearch_U<Real>(parlist) {
    const Real half(0.5);
    rho_ = parlist.sublist("Step").sublist("Line Search").sublist("Line-Search Method").get("Backtracking Rate",half);
  }

  void initialize(const Vector<Real> &x, const Vector<Real> &g) override {
    LineSearch_U<Real>::initialize(x,g);
    xnew_ = x.clone();
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
    // Perform backtracking
    while ( !status(LINESEARCH_U_BACKTRACKING,ls_neval,ls_ngrad,alpha,fold,gs,fval,*xnew_,s,obj) ) {
      alpha *= rho_;
      // Update iterate
      xnew_->set(x); xnew_->axpy(alpha,s);
      // Get objective value at xnew
      obj.update(*xnew_,UpdateType::Trial);
      fval = obj.value(*xnew_,tol);
      ls_neval++;
    }
  }
}; // class ROL::BackTracking_U

} // namespace ROL

#endif
