// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STATUSTEST_H
#define ROL_STATUSTEST_H

#include "ROL_Types.hpp"
#include "ROL_ParameterList.hpp"

/** \class ROL::StatusTest
    \brief Provides an interface to check status of optimization algorithms.
*/


namespace ROL {

template <class Real>
class StatusTest {
private:

  Real gtol_;
  Real stol_;
  int  max_iter_;
  bool use_rel_;

public:

  virtual ~StatusTest() {}

  StatusTest( ROL::ParameterList &parlist ) {
    Real em6(1e-6);
    gtol_     = parlist.sublist("Status Test").get("Gradient Tolerance", em6);
    stol_     = parlist.sublist("Status Test").get("Step Tolerance", em6*gtol_);
    max_iter_ = parlist.sublist("Status Test").get("Iteration Limit", 100);
    use_rel_  = parlist.sublist("Status Test").get("Use Relative Tolerances", false);
  }

  StatusTest( Real gtol = 1.e-6, Real stol = 1.e-12, int max_iter = 100, bool use_rel = false ) :  
    gtol_(gtol), stol_(stol), max_iter_(max_iter), use_rel_(use_rel) {}

  /** \brief Check algorithm status.
  */
  virtual bool check( AlgorithmState<Real> &state ) {
    if (state.iter==0 && use_rel_) {
      gtol_ *= state.gnorm;
      stol_ *= state.gnorm;
    }
    if ( (state.gnorm > gtol_) &&
         (state.snorm > stol_) &&
         (state.iter  < max_iter_) ) {
      return true;
    }
    else {
      state.statusFlag = (state.gnorm <= gtol_ ? EXITSTATUS_CONVERGED
                          : state.snorm <= stol_ ? EXITSTATUS_STEPTOL
                          : state.iter >= max_iter_ ? EXITSTATUS_MAXITER
                          : std::isnan(state.gnorm)||std::isnan(state.snorm) ? EXITSTATUS_NAN
                          : EXITSTATUS_LAST);
      return false;
    }
  }

}; // class StatusTest

} // namespace ROL

#endif
