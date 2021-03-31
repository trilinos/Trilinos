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
#ifndef ROL2_TYPEU_GRADIENT_H
#define ROL2_TYPEU_GRADIENT_H

#include "ROL2_TypeU_DescentDirection.hpp"
#include "ROL_Types.hpp"
#include "ROL_Step.hpp"

/** @ingroup step_group
    \class ROL::ROL2::TypeU::Gradient
    \brief Provides the interface to compute optimization steps
           with the gradient descent method globalized using line search.
*/

namespace ROL {
namespace ROL2 {
namespace TypeU {

template<typename Real>
class Gradient : public DescentDirection<Real> {
public:
  Gradient() {}

  void compute( Vector<Real> &s, Real &snorm, Real &sdotg, int &iter, int &flag,
          const Vector<Real> &x, const Vector<Real> &g, Objective<Real> &obj) override {
    s.set(g.dual());
    s.scale(static_cast<Real>(-1));
    snorm = s.norm();
    //sdotg = s.dot(g.dual());
    sdotg = s.apply(g);
    iter  = 0;
    flag  = 0;
  }

  void writeName( std::stream& os ) const override {
    os <<"Gradient Descent";
  }
}; // class ROL::ROL2::TypeU::Gradient
} // namespace TypeU 
} // namespace ROL2 
} // namespace ROL

namespace ROL2 {
namespace TypeU {
using ROL::ROL2::TypeU::Gradient
} // namespace TypeU 
} // namespace ROL2
#endif
