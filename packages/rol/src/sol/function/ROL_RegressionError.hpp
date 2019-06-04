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
