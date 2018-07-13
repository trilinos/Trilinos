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

#ifndef ROL_COMBINEDSTATUSTEST_H
#define ROL_COMBINEDSTATUSTEST_H

#include "ROL_StatusTest.hpp"

/** \class ROL::CombinedStatusTest
    \brief Provides an interface to check two status tests of optimization algorithms.
*/


namespace ROL {


template <class Real>
class CombinedStatusTest : public StatusTest<Real> {
private:
  std::vector<ROL::Ptr<StatusTest<Real> > > status_;

public:
  CombinedStatusTest(void) {
    status_.clear();
  }

  void reset(void) {
    status_.clear();
  }

  void add(const ROL::Ptr<StatusTest<Real> > &status) {
    status_.push_back(status);
  }

  bool check( AlgorithmState<Real> &state ) {
    int size = static_cast<int>(status_.size());
    if (size==0) {
      throw Exception::NotImplemented(">>> ROL::CombinedStatusTest::check: No status test has been added!");
    }

    bool flag = true;
    for (int i = 0; i < size; ++i) {
      if (!(status_[i]->check(state))) {
        flag = false;
        break;
      }
    }

    // true  = "not converged"
    // false = "converged"
    return flag;
  }

}; // class CombinedStatusTest

} // namespace ROL

#endif
