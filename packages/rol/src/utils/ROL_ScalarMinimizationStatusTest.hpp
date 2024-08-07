// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SCALARMINIMIZATIONSTATUSTEST_H
#define ROL_SCALARMINIMIZATIONSTATUSTEST_H

/** \class ROL::ScalarMinimizationStatusTest
    \brief Provides interface for terminating scalar minimization algoithms.
*/

namespace ROL { 

template<class Real>
class ScalarMinimizationStatusTest {
public:
  virtual ~ScalarMinimizationStatusTest() {}
  virtual bool check(Real &x, Real &fx, Real &gx,
                     int &nfval, int &ngrad,
                     const bool deriv = false) {
    return false;
  }
};

}

#endif
