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
    const std::vector<Real> data = Objective<Real>::getParameter();
    const unsigned dim = x.size();
    Real val = data[0] - x[0];
    for (unsigned i = 1; i < dim; ++i) val -= data[i] * x[i];
    return val;
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    checkSize(g);
    const std::vector<Real> data = Objective<Real>::getParameter();
    const unsigned dim = g.size();
    g[0] = static_cast<Real>(-1);
    for (unsigned i = 1; i < dim; ++i) g[i] = -data[i];
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    checkSize(hv);
    hv.assign(hv.size(),static_cast<Real>(0));
  }
}; // class RegressionError

} // namespace ROL

#endif
