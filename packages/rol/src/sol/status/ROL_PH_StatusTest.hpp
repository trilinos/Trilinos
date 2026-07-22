// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PH_STATUSTEST_H
#define ROL_PH_STATUSTEST_H

#include "ROL_StatusTest.hpp"

/** \class ROL::PH_StatusTest
    \brief Provides an interface to check status of the progressive hedging algorithm.
*/


namespace ROL {

template <class Real>
class PH_StatusTest : public StatusTest<Real> {
private:

  Real mu_;
  Real epsilon_;
  Ptr<const Vector<Real>> xbar_;
  Real tol_;
  Ptr<Vector<Real>> x_;

public:

  PH_StatusTest( ROL::ParameterList &parlist, const Vector<Real> &x ) {
    mu_      = parlist.sublist("SOL").sublist("Progressive Hedging").get("Fixed Tolerance", 1e-4);
    epsilon_ = parlist.sublist("SOL").sublist("Progressive Hedging").get("Dynamic Tolerance", 0.1);
    x_       = x.clone();
  }

  void setData(const int iter, const Ptr<const Vector<Real>> &xbar) {
    const Real one(1);
    tol_   = mu_*std::pow(one-epsilon_,iter+1);
    xbar_  = xbar;
  }

  bool check( AlgorithmState<Real> &state ) {
    const Real one(1);
    x_->set(*state.iterateVec); x_->axpy(-one,*xbar_);
    Real xnorm = x_->norm();
    if ( state.gnorm <= tol_*std::min(one,xnorm) ) {
      state.statusFlag = EXITSTATUS_USERDEFINED;
      return false;
    }
    return true;
  }

}; // class PH_StatusTest

} // namespace ROL

#endif
