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

/*! \file  pde_ginzburg-landau.hpp
    \brief Implements the local PDE interface for the optimal control of
           simplified Ginzburg-Landau.
*/

#ifndef PDE_GINZBURGLANDAU_EX02_HPP
#define PDE_GINZBURGLANDAU_EX02_HPP

#include "pde_ginzburg-landau.hpp"


template <class Real>
class PDE_GinzburgLandau_ex02 : public PDE_GinzburgLandau<Real> {
public:
  PDE_GinzburgLandau_ex02(Teuchos::ParameterList &parlist) : PDE_GinzburgLandau<Real>(parlist) {}

  void evaluateMagneticPotential(std::vector<Real> &Ax, const std::vector<Real> &x) const {
    const Real pi(M_PI);
    Ax[0] =  std::sin(pi*x[0])*std::cos(pi*x[1]);
    Ax[1] = -std::cos(pi*x[0])*std::sin(pi*x[1]);
  }

  Real evaluateNeumann(const std::vector<Real> &x, const int component) const {
    return static_cast<Real>(0);
  }

  Real evaluateForce(const std::vector<Real> &x, const int component) const {
    const Real pi(M_PI), one(1), two(2);
    const Real cx = std::cos(pi*x[0]), sx = std::sin(pi*x[0]);
    const Real cy = std::cos(pi*x[1]), sy = std::sin(pi*x[1]);
    const Real pi2 = pi*pi, cx2 = cx*cx, cy2 = cy*cy, sx2 = sx*sx, sy2 = sy*sy;
    return (component == 0) ? pi2*cx + (sx2*cy2 + cx2*sy2 - one)*cx
                              -two*pi*cx*sy2 + (cx2 + cy2)*cx
                            : pi2*cy + (sx2*cy2 + cx2*sy2 - one)*cy
                              -two*pi*sx2*cy + (cx2 + cy2)*cy;
  }

}; // PDE_GinzburgLandau_ex01

#endif
