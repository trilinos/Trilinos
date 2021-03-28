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

#pragma once
#ifndef ROL2_STATUSTEST_DEF_H
#define ROL2_STATUSTEST_DEF_H

/** \class ROL2::StatusTest
    \brief Provides an interface to check status of optimization algorithms.
*/


namespace ROL2 {

template<class Real>
StatusTest<Real>::StatusTest( ParameterList& parlist ) {
  Real em6(1e-6);
  auto& stlist = parlist.sublist("Status Test");
  gtol_     = stlist.get("Gradient Tolerance", em6);
  stol_     = stlist.get("Step Tolerance", em6*gtol_);
  max_iter_ = stlist.get("Iteration Limit", 100);
  use_rel_  = stlist.get("Use Relative Tolerances", false);
}

template<class Real>
StatusTest<Real>::StatusTest( Real gtol, 
                              Real stol, 
                              int  max_iter, 
                              bool use_rel ) 
  : gtol_(gtol), stol_(stol), max_iter_(max_iter), use_rel_(use_rel) {}

template<class Real>
bool StatusTest<Real>::check( typename Algorithm<Real>::State& state ) {

   if (state.iter_ ==0 && use_rel_) {
     gtol_ *= state.gnorm_;
     stol_ *= state.gnorm_;
   }
   if ( (state.gnorm_ > gtol_) && 
        (state.snorm_ > stol_) && 
        (state.iter_  < max_iter_) ) {
     return true;
   }
   else {
     state.statusFlag_ = (state.gnorm_ <= gtol_ ? ExitStatus::Converged
                         : state.snorm_ <= stol_ ? ExitStatus::StepTolMet
                         : state.iter_ >= max_iter_ ? ExitStatus::MaxIter
                         : std::isnan(state.gnorm_)||std::isnan(state.snorm_) ? ExitStatus::EncounteredNaN
                         : ExitStatus::Last);
     return false;
   }
} // StatusTest<Real>::check

} // namespace ROL2

#endif // ROL2_STATUSTEST_DEF_HPP
