// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CONSTRAINTSTATUSTEST_H
#define ROL_CONSTRAINTSTATUSTEST_H

#include "ROL_StatusTest.hpp"

/** \class ROL::ConstraintStatusTest
    \brief Provides an interface to check status of optimization algorithms
           for problems with equality constraints.
*/


namespace ROL {

template <class Real>
class ConstraintStatusTest : public StatusTest<Real> {
private:

  Real gtol_, gtol0_;
  Real ctol_, ctol0_;
  Real stol_, stol0_;
  int  max_iter_;
  bool use_rel_;

public:

  virtual ~ConstraintStatusTest() {}

  ConstraintStatusTest( ROL::ParameterList &parlist ) {
    Real em6(1e-6);
    gtol_     = parlist.sublist("Status Test").get("Gradient Tolerance", em6);
    ctol_     = parlist.sublist("Status Test").get("Constraint Tolerance", em6);
    stol_     = parlist.sublist("Status Test").get("Step Tolerance", em6*gtol_);
    max_iter_ = parlist.sublist("Status Test").get("Iteration Limit", 100);
    use_rel_  = parlist.sublist("Status Test").get("Use Relative Tolerances", false);
    gtol0_    = gtol_;
    ctol0_    = ctol_;
    stol0_    = stol_;
  }

  ConstraintStatusTest( Real gtol = 1e-6, Real ctol = 1e-6, Real stol = 1e-12, int max_iter = 100, bool use_rel = false ) :
    gtol_(gtol), gtol0_(gtol), ctol_(ctol), ctol0_(ctol), stol_(stol), stol0_(stol), max_iter_(max_iter), use_rel_(use_rel) {}

  /** \brief Check algorithm status.
  */
  virtual bool check( AlgorithmState<Real> &state ) {
    if (state.iter==0 && use_rel_) {
      gtol_ = gtol0_*std::max(state.gnorm,static_cast<Real>(1e-2));
      ctol_ = ctol0_*std::max(state.cnorm,static_cast<Real>(1e-2));
      stol_ = stol0_*std::max(std::min(state.gnorm,state.cnorm),static_cast<Real>(1e-2));
    }
    if ( ((state.gnorm > gtol_) || (state.cnorm > ctol_)) && 
          (state.snorm > stol_) && 
          (state.iter  < max_iter_) ) {
      return true;
    }
    else {
      state.statusFlag = ((state.gnorm <= gtol_) && (state.cnorm <= ctol_) ? EXITSTATUS_CONVERGED
                           : state.snorm <= stol_ ? EXITSTATUS_STEPTOL
                           : state.iter  >= max_iter_ ? EXITSTATUS_MAXITER
                           : EXITSTATUS_LAST);
      return false;
    }
  }

}; // class ConstraintStatusTest

} // namespace ROL

#endif
