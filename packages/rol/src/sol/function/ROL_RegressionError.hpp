// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_REGRESSIONERROR_H
#define ROL_REGRESSIONERROR_H

#include "ROL_StdVector.hpp"
#include "ROL_StdObjective.hpp"

/** @ingroup func_group
    \class ROL::RegressionError
    \brief Provides the interface to evaluate linear regression error.

    ---
*/


namespace ROL {

template <class Real>
class RegressionError : public StdObjective<Real> {
private:
  void checkSize(const std::vector<Real> &x) {
    const std::vector<Real> data = Objective<Real>::getParameter();
    if (data.size() != x.size()) {
      throw Exception::NotImplemented("ROL::RegressionError : Data dimension does not match input dimension!");
    }
  }

public:
  RegressionError(void) {}

  Real value( const std::vector<Real> &x, Real &tol ) {
    checkSize(x);
    // Parse Input Vector
    std::vector<Real> c; c.assign(x.begin()+1,x.end());
    Real c0 = x[0];
    // Parse Data
    const std::vector<Real> data = Objective<Real>::getParameter();
    std::vector<Real> X; X.assign(data.begin()+1,data.end());
    Real Y = data[0];
    // Build Error
    int Xdim = X.size();
    Real val = Y-c0;
    for (int i = 0; i < Xdim; ++i) {
      val -= c[i]*X[i];
    }
    return val;
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    checkSize(x); checkSize(g);
    // Parse Data
    const std::vector<Real> data = Objective<Real>::getParameter();
    std::vector<Real> X; X.assign(data.begin()+1,data.end());
    // Build Error
    int Xdim = X.size();
    g[0] = static_cast<Real>(-1);
    for (int i = 0; i < Xdim; ++i) {
      g[i+1] = -X[i];
    }
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    hv.assign(hv.size(),static_cast<Real>(0));
  }
}; // class RegressionError

} // namespace ROL

#endif
