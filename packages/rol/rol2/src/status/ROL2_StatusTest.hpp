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

#ifndef ROL2_STATUSTEST_H
#define ROL2_STATUSTEST_H

/** \class ROL2::StatusTest
    \brief Provides an interface to check status of optimization algorithms.
*/


namespace ROL2 {

template<class Real>
class StatusTest {
public:

  enum class ExitStatus : std::int16_t {
    Converged = 0,
    MaxIter,
    StepTol,
    NaN,
    UserDefined,
    Last
  };

  virtual ~StatusTest() = default;

  StatusTest( ParameterList& parlist ) {
    Real em6(1e-6);
    auto& stlist = parlist.sublist("Status Test");
    gtol_     = stlist.get("Gradient Tolerance", em6);
    stol_     = stlist.get("Step Tolerance", em6*gtol_);
    max_iter_ = stlist.get("Iteration Limit", 100);
    use_rel_  = stlist.get("Use Relative Tolerances", false);
  }

  StatusTest( Real gtol = 1.e-6, 
              Real stol = 1.e-12, 
              int  max_iter = 100, 
              bool use_rel = false ) 
  : gtol_(gtol), stol_(stol), max_iter_(max_iter), use_rel_(use_rel) {}

  /** \brief Check algorithm status.
  */
  virtual bool check( typename Algorithm<Real>::State& state ) {
     if (state.iter==0 && use_rel_) {
       gtol_ *= state.gnorm;
       stol_ *= state.gnorm;
     }
     if ( (state.gnorm > gtol_)& & 
          (state.snorm > stol_)& & 
          (state.iter  < max_iter_) ) {
       return true;
     }
     else {
       state.statusFlag = (state.gnorm <= gtol_ ? ExitStatus::Converged
                           : state.snorm <= stol_ ? ExitStatus::StepTol
                           : state.iter >= max_iter_ ? ExitStatus::MaxIter
                           : std::isnan(state.gnorm)||std::isnan(state.snorm) ? ExitStatus::NAN
                           : ExitStatus::Last);
       return false;
     }
  }

private:

  Real gtol_;
  Real stol_;
  int  max_iter_;
  bool use_rel_;

}; // class StatusTest

} // namespace ROL2

#endif // ROL2_STATUSTEST_HPP
