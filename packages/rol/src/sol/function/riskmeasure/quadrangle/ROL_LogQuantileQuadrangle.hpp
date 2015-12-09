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

#ifndef ROL_LOGQUANTILEQUAD_HPP
#define ROL_LOGQUANTILEQUAD_HPP

#include "ROL_ExpectationQuad.hpp"
#include "ROL_PlusFunction.hpp"

namespace ROL {

template<class Real>
class LogQuantileQuadrangle : public ExpectationQuad<Real> {
private:

  Teuchos::RCP<PlusFunction<Real> > pf_;

  Real prob_;
  Real eps_;

public:

  LogQuantileQuadrangle(Real prob, Real eps, Teuchos::RCP<PlusFunction<Real> > &pf ) 
    : ExpectationQuad<Real>(), pf_(pf) {
    prob_ = ((prob >= 0.0) ? ((prob <= 1.0) ? prob : 0.5) : 0.5);
    eps_  = ((eps > 0.0) ? eps : 1.0);
  }

  LogQuantileQuadrangle(Teuchos::ParameterList &parlist) : ExpectationQuad<Real>() {
    Teuchos::ParameterList& list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Log-Quantile Quadrangle");
    // Check CVaR inputs
    Real prob = list.get("Confidence Level",0.5);
    prob_ = ((prob >= 0.0) ? ((prob <= 1.0) ? prob : 0.5) : 0.5);
    // Build plus function
    pf_ = Teuchos::rcp( new PlusFunction<Real>(list) );
    Real eps = list.get("Smoothing Parameter",1.);
    eps_ = ((eps > 0.) ? eps : 1.);
  }

  Real error(Real x, int deriv = 0) {
    TEUCHOS_TEST_FOR_EXCEPTION( (deriv > 2), std::invalid_argument,
      ">>> ERROR (ROL::LogQuantileQuadrangle::error): deriv greater than 2!");
    TEUCHOS_TEST_FOR_EXCEPTION( (deriv < 0), std::invalid_argument,
      ">>> ERROR (ROL::LogQuantileQuadrangle::error): deriv less than 0!");

    Real X = ((deriv == 0) ? x : ((deriv == 1) ? 1.0 : 0.0));
    return regret(x,deriv) - X;
  }

  Real regret(Real x, int deriv = 0) {
    TEUCHOS_TEST_FOR_EXCEPTION( (deriv > 2), std::invalid_argument,
      ">>> ERROR (ROL::LogQuantileQuadrangle::regret): deriv greater than 2!");
    TEUCHOS_TEST_FOR_EXCEPTION( (deriv < 0), std::invalid_argument,
      ">>> ERROR (ROL::LogQuantileQuadrangle::regret): deriv less than 0!");

    Real arg = std::exp(x);
    Real reg = 1.0/(1.0-prob_) * (pf_->evaluate(arg-1.0,deriv) *
                 ((deriv == 0) ? 1.0 : ((deriv == 1) ? arg : arg*arg))
               + ((deriv == 2) ? pf_->evaluate(arg-1.0,deriv-1)*arg : 0.0));
    return reg;
  }

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
