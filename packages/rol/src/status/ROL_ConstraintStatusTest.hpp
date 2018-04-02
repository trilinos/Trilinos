// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
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

  Real gtol_;
  Real ctol_;
  Real stol_;
  int  max_iter_;

public:

  virtual ~ConstraintStatusTest() {}

  ConstraintStatusTest( Teuchos::ParameterList &parlist ) {
    Real em6(1e-6);
    gtol_     = parlist.sublist("Status Test").get("Gradient Tolerance", em6);
    ctol_     = parlist.sublist("Status Test").get("Constraint Tolerance", em6);
    stol_     = parlist.sublist("Status Test").get("Step Tolerance", em6*gtol_);
    max_iter_ = parlist.sublist("Status Test").get("Iteration Limit", 100);
  }

  ConstraintStatusTest( Real gtol = 1e-6, Real ctol = 1e-6, Real stol = 1e-12, int max_iter = 100 ) :  
    gtol_(gtol), ctol_(ctol), stol_(stol), max_iter_(max_iter) {}

  /** \brief Check algorithm status.
  */
  virtual bool check( AlgorithmState<Real> &state ) {
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
