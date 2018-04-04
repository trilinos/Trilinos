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

  BundleStatusTest( Teuchos::ParameterList &parlist ) {
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
