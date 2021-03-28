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

#ifndef ROL2_COMBINEDSTATUSTEST_HPP
#define ROL2_COMBINEDSTATUSTEST_HPP

/** \class ROL2::CombinedStatusTest
    \brief Provides an interface to check two status tests of optimization algorithms.
*/

namespace ROL2 {

template<typename Real>
CombinedStatusTest<Real>::CombinedStatusTest() {
  status_.clear();
}

template<typename Real>
void CombinedStatusTest<Real>::reset() {
  status_.clear();
}

template<typename Real>
void CombinedStatusTest<Real>::add( const Ptr<StatusTest<Real>>& status ) {
  status_.push_back(status);
}

template<typename Real>
bool CombinedStatusTest<Real>::check( typename Algorithm<Real>::State& state ) {
  ROL_TEST_FOR_EXCEPTION(status_.empty(),std::logic_error,
    ">>> ROL::CombinedStatusTest::check : No status test has been added!");

  bool flag = true;

  for ( const auto& status : status_) {
    flag = status->check(state);
    if (!flag) break;
  }

  return flag;
}

} // namespace ROL2

#endif // ROL2_COMBINEDSTATUSTEST_HPP
