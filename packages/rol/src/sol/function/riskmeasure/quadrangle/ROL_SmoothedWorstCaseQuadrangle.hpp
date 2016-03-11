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

#ifndef ROL_SMOOTHEDWORSTCASEQUAD_HPP
#define ROL_SMOOTHEDWORSTCASEQUAD_HPP

#include "ROL_ExpectationQuad.hpp"

namespace ROL {

template<class Real>
class SmoothedWorstCaseQuadrangle : public ExpectationQuad<Real> {
private:

  Real eps_;

public:

  SmoothedWorstCaseQuadrangle(const Real eps) : ExpectationQuad<Real>() {
    Real zero(0), one(1);
    eps_  = ((eps > zero) ? eps : one);
  }

  SmoothedWorstCaseQuadrangle(Teuchos::ParameterList &parlist) : ExpectationQuad<Real>() {
    Real zero(0), one(1);
    Teuchos::ParameterList& list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Smoothed Worst-Case Quadrangle");
    Real eps = list.get("Smoothing Parameter",one);
    eps_  = ((eps > zero) ? eps : one);
  }

  Real error(Real x, int deriv = 0) {
    Real err(0), zero(0), half(0.5), one(1);
    if (deriv == 0) {
      err = (x <= -eps_) ? -(x+half*eps_) : half*x*x/eps_;
    }
    else if (deriv == 1) {
      err = (x <= -eps_) ? -one : x/eps_;
    }
    else {
      err = (x <= -eps_) ? zero : one/eps_;
    }
    return err;
  }

  Real regret(Real x, int deriv = 0) {
    Real zero(0), one(1);
    Real X = ((deriv==0) ? x : ((deriv==1) ? one : zero));
    Real reg = error(x,deriv) + X;
    return reg;
  }

  void checkRegret(void) {
    ExpectationQuad<Real>::checkRegret();
    // Check v'(-eps)
    Real x = -eps_, zero(0), one(1), two(2), p1(0.1);
    Real vx = zero, vy = zero;
    Real dv = regret(x,1);
    Real t = one;
    Real diff = zero;
    Real err = zero;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v'(-eps) is correct? \n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "v'(x)"
                            << std::setw(20) << "(v(x+t)-v(x-t))/2t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      vy = regret(x+t,0);
      vx = regret(x-t,0);
      diff = (vy-vx)/(two*t);
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << t
                << std::setw(20) << dv
                << std::setw(20) << diff
                << std::setw(20) << err
                << "\n";
      t *= p1;
    }
    std::cout << "\n";
    // check v''(-eps) 
    vx = zero;
    vy = zero;
    dv = regret(x,2);
    t = one;
    diff = zero;
    err = zero;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v''(-eps) is correct? \n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "v''(x)"
                            << std::setw(20) << "(v'(x+t)-v'(x-t))/2t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      vy = regret(x+t,1);
      vx = regret(x-t,1);
      diff = (vy-vx)/(two*t);
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << t
                << std::setw(20) << dv
                << std::setw(20) << diff
                << std::setw(20) << err
                << "\n";
      t *= p1;
    }
    std::cout << "\n"; 
    // Check v'(0)
    x = zero;
    vx = zero;
    vy = zero;
    dv = regret(x,1);
    t = one;
    diff = zero;
    err = zero;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v'(0) is correct? \n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "v'(x)"
                            << std::setw(20) << "(v(x+t)-v(x-t))/2t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      vy = regret(x+t,0);
      vx = regret(x-t,0);
      diff = (vy-vx)/(two*t);
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << t
                << std::setw(20) << dv
                << std::setw(20) << diff
                << std::setw(20) << err
                << "\n";
      t *= p1;
    }
    std::cout << "\n";
    // check v''(0) 
    vx = zero;
    vy = zero;
    dv = regret(x,2);
    t = one;
    diff = zero;
    err = zero;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v''(0) is correct? \n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "v''(x)"
                            << std::setw(20) << "(v'(x+t)-v'(x-t))/2t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      vy = regret(x+t,1);
      vx = regret(x-t,1);
      diff = (vy-vx)/(two*t);
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << t
                << std::setw(20) << dv
                << std::setw(20) << diff
                << std::setw(20) << err
                << "\n";
      t *= p1;
    }
    std::cout << "\n"; 
    // Check v'(-1-eps)
    x = -eps_-one;
    vx = zero;
    vy = zero;
    dv = regret(x,1);
    t = one;
    diff = zero;
    err = zero;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v'(-1-eps) is correct? \n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "v'(x)"
                            << std::setw(20) << "(v(x+t)-v(x-t))/2t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      vy = regret(x+t,0);
      vx = regret(x-t,0);
      diff = (vy-vx)/(two*t);
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << t
                << std::setw(20) << dv
                << std::setw(20) << diff
                << std::setw(20) << err
                << "\n";
      t *= p1;
    }
    std::cout << "\n";
    // check v''(-1-eps) 
    vx = zero;
    vy = zero;
    dv = regret(x,2);
    t = one;
    diff = zero;
    err = zero;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v''(-1-eps) is correct? \n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "v''(x)"
                            << std::setw(20) << "(v'(x+t)-v'(x-t))/2t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      vy = regret(x+t,1);
      vx = regret(x-t,1);
      diff = (vy-vx)/(two*t);
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << t
                << std::setw(20) << dv
                << std::setw(20) << diff
                << std::setw(20) << err
                << "\n";
      t *= p1;
    }
    std::cout << "\n"; 
  }

};

}
#endif
