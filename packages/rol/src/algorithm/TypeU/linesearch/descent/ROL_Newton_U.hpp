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

#ifndef ROL_NEWTON_U_H
#define ROL_NEWTON_U_H

#include "ROL_DescentDirection_U.hpp"
#include "ROL_Types.hpp"

/** @ingroup step_group
    \class ROL::Newton_U
    \brief Provides the interface to compute optimization steps
           with Newton's method globalized using line search.
*/

namespace ROL {

template<typename Real>
class Newton_U : public DescentDirection_U<Real> {
public:
  Newton_U() {}

  void compute( Vector<Real> &s, Real &snorm, Real &sdotg, int &iter, int &flag,
          const Vector<Real> &x, const Vector<Real> &g, Objective<Real> &obj) override {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    // Compute unconstrained step
    obj.invHessVec(s,g,x,tol);
    //sdotg = -s.dot(g.dual());
    sdotg = -s.apply(g);
    if (sdotg >= static_cast<Real>(0)) {
      s.set(g.dual());
      //sdotg = -s.dot(g.dual());
      sdotg = -s.apply(g);
    }
    s.scale(static_cast<Real>(-1));
    snorm = s.norm();
    iter  = 0;
    flag  = 0;
  }

  std::string printName(void) const override {
    std::string name = "Newton's Method";
    return name;
  }
}; // class ROL::TypeU::LineSearch::Newton

} // namespace ROL

#endif
