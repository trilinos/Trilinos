// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BUNDLESTATUSTEST_H
#define ROL_BUNDLESTATUSTEST_H

#include "ROL_StatusTest.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template <class Real>
class BundleStatusTest : public StatusTest<Real> {
private:

  Real tol_;
  int  max_iter_;

public:

  virtual ~BundleStatusTest() {}

  BundleStatusTest( ROL::ParameterList &parlist ) {
    Real em6(1e-6);
    tol_      = parlist.sublist("Step").sublist("Bundle").get("Epsilon Solution Tolerance", em6);
    max_iter_ = parlist.sublist("Status Test").get("Iteration Limit", 100);
  }

  BundleStatusTest( Real tol = 1.e-6, int max_iter = 100 ) :
    tol_(tol), max_iter_(max_iter) {}

  /** \brief Check algorithm status.
  */
  virtual bool check( AlgorithmState<Real> &state ) {
     bool stat = false;
     if ( (std::max(state.aggregateGradientNorm,state.aggregateModelError) > tol_)  
         && (state.iter < max_iter_) 
         && (state.flag == false) ) {
       stat = true;
     }
     else {
       state.statusFlag = (std::max(state.aggregateGradientNorm,state.aggregateModelError) <= tol_ ? EXITSTATUS_CONVERGED
                           : state.iter >= max_iter_ ? EXITSTATUS_MAXITER
                           : state.flag == true ? EXITSTATUS_CONVERGED
                           : EXITSTATUS_LAST);
     }
     return stat;
  }

}; // class BundleStatusTest

} // namespace ROL

#endif
