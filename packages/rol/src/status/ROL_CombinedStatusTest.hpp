// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_COMBINEDSTATUSTEST_H
#define ROL_COMBINEDSTATUSTEST_H

#include "ROL_StatusTest.hpp"

/** \class ROL::CombinedStatusTest
    \brief Provides an interface to check two status tests of optimization algorithms.
*/


namespace ROL {

template<typename Real>
class CombinedStatusTest : public StatusTest<Real> {
private:
  std::vector<Ptr<StatusTest<Real>>> status_;

public:
  CombinedStatusTest(void) {
    status_.clear();
  }

  void reset(void) {
    status_.clear();
  }

  void add(const Ptr<StatusTest<Real>> &status) {
    status_.push_back(status);
  }

  bool check( AlgorithmState<Real> &state ) {
    ROL_TEST_FOR_EXCEPTION(status_.empty(),std::logic_error,
      ">>> ROL::CombinedStatusTest::check : No status test has been added!");

    bool flag = true;
    for (const auto & status : status_) {
      flag = status->check(state);
      if (!flag) break;
    }

    // true  = "continue iteration"
    // false = "stop iteration"
    return flag;
  }

}; // class CombinedStatusTest

} // namespace ROL

#endif
