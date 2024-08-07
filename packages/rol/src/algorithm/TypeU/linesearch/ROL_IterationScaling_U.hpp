// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_ITERATIONSCALING_U_H
#define ROL_ITERATIONSCALING_U_H

/** \class ROL::IterationScaling_U
    \brief Provides an implementation of iteration scaled line search.
*/

#include "ROL_LineSearch_U.hpp"

namespace ROL {

template<typename Real>
class IterationScaling_U : public LineSearch_U<Real> {
private:
  int algo_iter_;
  Ptr<Vector<Real>> xnew_; 

  using LineSearch_U<Real>::getInitialAlpha;

public:

  IterationScaling_U( ParameterList &parlist ) : LineSearch_U<Real>(parlist), algo_iter_(0) {}

  void initialize(const Vector<Real> &x, const Vector<Real> &g) override {
    LineSearch_U<Real>::initialize(x,g);
    xnew_   = x.clone();
  }

  // Run Iteration scaled line search
  void run( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
            const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
            Objective<Real> &obj ) override {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    ls_neval = 0;
    ls_ngrad = 0;
    // Get line search parameter
    algo_iter_++;
    alpha = getInitialAlpha(ls_neval,ls_ngrad,fval,gs,x,s,obj)/static_cast<Real>(algo_iter_);
    // Update iterate
    xnew_->set(x); xnew_->axpy(alpha,s);
    // Compute objective function value
    obj.update(*xnew_,UpdateType::Trial);
    fval = obj.value(*xnew_,tol);
    ls_neval++;
  }
}; // class ROL::IterationScaling_U

} // namespace ROL

#endif
