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

#ifndef ROL_SMOOTHCVARQUAD_HPP
#define ROL_SMOOTHCVARQUAD_HPP

#include "ROL_ExpectationQuad.hpp"
#include "ROL_PlusFunction.hpp"

namespace ROL {

template<class Real>
class SmoothCVaRQuad : public ExpectationQuad<Real> {
private:
  Real prob_;
  Real eps_;
  Teuchos::RCP<PlusFunction<Real> > pf_;

public:
  SmoothCVaRQuad(Real prob, Real eps, Teuchos::RCP<PlusFunction<Real> > &pf ) 
    : prob_(prob), eps_(eps), pf_(pf), ExpectationQuad<Real>() {}

  Real error(Real x, int deriv = 0) {
    Real err = (prob_/(1.0-prob_))*pf_->evaluate(x,deriv) 
              + ((deriv%2) ? -1.0 : 1.0)*pf_->evaluate(-x,deriv);
    return err;
  }

  Real regret(Real x, int deriv = 0) {
    Real X = ((deriv==0) ? x : ((deriv==1) ? 1.0 : 0.0));
    Real reg = error(x,deriv) + X;
    return reg;
  }

/*
  // This is derived from a smoothed Koenker-Bassett error function
  Real regret(Real x, int deriv = 0) {
    int dom = ((x >= eps_) ? 0 : ((x < eps_ && x >= 0.0) ? 1 : ((x < 0.0 && x > -eps_) ? 2 : 3)));

    Real val = 0.0;
    Real c1  = prob_/(1.0-prob_);
    Real c2  = 1.0/(1.0-prob_);
    Real x2  = x*x;
    Real x3  = x*x2;
    Real x4  = x*x3;
    Real e1  = eps_*0.5;
    Real e2  = eps_*eps_;
    Real e3  = e2*eps_*2.0;
    switch (dom) {
      case 0: // x is greater than or equal to eps 
        switch (deriv) {
          case 0: val = c2*x - c1*e1; break;
          case 1: val = c2;           break;
          case 2: val = 0.0;          break;
        }
        break;
      case 1: // x is between 0 and eps
        switch (deriv) {
          case 0: val = c1*(x3/e2 - x4/e3) + x;           break;
          case 1: val = c1*(3.0*x2/e2 - 4.0*x3/e3) + 1.0; break;
          case 2: val = c1*(6.0*x/e2 - 12.0*x2/e3);       break;
        }
        break;
      case 2: // x is between -eps and 0
        switch (deriv) {
          case 0: val = (-x3/e2 - x4/e3) + x;           break;
          case 1: val = (-3.0*x2/e2 - 4.0*x3/e3) + 1.0; break;
          case 2: val = (-6.0*x/e2 - 12.0*x2/e3);       break;
        }
        break;
      case 3: // x is less than or equal to eps
        switch (deriv) {
          case 0: val = -e1; break;
          case 1: val = 0.0; break;
          case 2: val = 0.0; break;
        }
        break;
    }
    return val;
  }
*/

  void checkRegret(void) {
    ExpectationQuad<Real>::checkRegret();
    // Check v'(eps)
    Real x = eps_;
    Real vx = 0.0, vy = 0.0;
    Real dv = regret(x,1);
    Real t = 1.0;
    Real diff = 0.0;
    Real err = 0.0;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v'(eps) is correct? \n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "v'(x)"
                            << std::setw(20) << "(v(x+t)-v(x-t))/2t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      vy = regret(x+t,0);
      vx = regret(x-t,0);
      diff = (vy-vx)/(2.0*t);
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << t
                << std::setw(20) << dv
                << std::setw(20) << diff
                << std::setw(20) << err
                << "\n";
      t *= 0.1;
    }
    std::cout << "\n";
    // check v''(eps) 
    vx = 0.0;
    vy = 0.0;
    dv = regret(x,2);
    t = 1.0;
    diff = 0.0;
    err = 0.0;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v''(eps) is correct? \n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "v''(x)"
                            << std::setw(20) << "(v'(x+t)-v'(x-t))/2t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      vy = regret(x+t,1);
      vx = regret(x-t,1);
      diff = (vy-vx)/(2.0*t);
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << t
                << std::setw(20) << dv
                << std::setw(20) << diff
                << std::setw(20) << err
                << "\n";
      t *= 0.1;
    }
    std::cout << "\n"; 
    // Check v'(0)
    x = 0.0;
    vx = 0.0;
    vy = 0.0;
    dv = regret(x,1);
    t = 1.0;
    diff = 0.0;
    err = 0.0;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v'(0) is correct? \n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "v'(x)"
                            << std::setw(20) << "(v(x+t)-v(x-t))/2t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      vy = regret(x+t,0);
      vx = regret(x-t,0);
      diff = (vy-vx)/(2.0*t);
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << t
                << std::setw(20) << dv
                << std::setw(20) << diff
                << std::setw(20) << err
                << "\n";
      t *= 0.1;
    }
    std::cout << "\n";
    // check v''(eps) 
    vx = 0.0;
    vy = 0.0;
    dv = regret(x,2);
    t = 1.0;
    diff = 0.0;
    err = 0.0;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v''(0) is correct? \n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "v''(x)"
                            << std::setw(20) << "(v'(x+t)-v'(x-t))/2t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      vy = regret(x+t,1);
      vx = regret(x-t,1);
      diff = (vy-vx)/(2.0*t);
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << t
                << std::setw(20) << dv
                << std::setw(20) << diff
                << std::setw(20) << err
                << "\n";
      t *= 0.1;
    }
    std::cout << "\n"; 
    // Check v'(0)
    x = -eps_;
    vx = 0.0;
    vy = 0.0;
    dv = regret(x,1);
    t = 1.0;
    diff = 0.0;
    err = 0.0;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v'(-eps) is correct? \n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "v'(x)"
                            << std::setw(20) << "(v(x+t)-v(x-t))/2t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      vy = regret(x+t,0);
      vx = regret(x-t,0);
      diff = (vy-vx)/(2.0*t);
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << t
                << std::setw(20) << dv
                << std::setw(20) << diff
                << std::setw(20) << err
                << "\n";
      t *= 0.1;
    }
    std::cout << "\n";
    // check v''(eps) 
    vx = 0.0;
    vy = 0.0;
    dv = regret(x,2);
    t = 1.0;
    diff = 0.0;
    err = 0.0;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v''(-eps) is correct? \n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "v''(x)"
                            << std::setw(20) << "(v'(x+t)-v'(x-t))/2t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      vy = regret(x+t,1);
      vx = regret(x-t,1);
      diff = (vy-vx)/(2.0*t);
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << t
                << std::setw(20) << dv
                << std::setw(20) << diff
                << std::setw(20) << err
                << "\n";
      t *= 0.1;
    }
    std::cout << "\n"; 
  }

};

}
#endif
