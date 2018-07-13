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

#ifndef ROL_MOREAUYOSIDACVAR_HPP
#define ROL_MOREAUYOSIDACVAR_HPP

#include "ROL_ExpectationQuad.hpp"

/** @ingroup risk_group
    \class ROL::MoreauYosidaCVaR
    \brief Provides an interface for a smooth approximation of the conditional
           value-at-risk.

    The conditional value-at-risk (also called the average value-at-risk
    or the expected shortfall) with confidence level \f$0\le \beta < 1\f$
    is
    \f[
       \mathcal{R}(X) = \inf_{t\in\mathbb{R}} \left\{
         t + \frac{1}{1-\beta} \mathbb{E}\left[(X-t)_+\right]
         \right\}
    \f]
    where \f$(x)_+ = \max\{0,x\}\f$.  If the distribution of \f$X\f$ is
    continuous, then \f$\mathcal{R}\f$ is the conditional expectation of
    \f$X\f$ exceeding the \f$\beta\f$-quantile of \f$X\f$ and the optimal
    \f$t\f$ is the \f$\beta\f$-quantile.
    Additionally, \f$\mathcal{R}\f$ is a law-invariant coherent risk measure.

    The conditional value-at-risk is in general not smooth due to
    \f$(\cdot)_+\f$.  One approach to smoothing the conditional-value-at-risk
    is to regularize its biconjugate form.  That is, since \f$\mathcal{R}\f$
    is coherent, we have that
    \f[
       \mathcal{R}(X) = \sup_{\vartheta\in\mathfrak{A}}\mathbb{E}[\vartheta X]
    \f]
    where \f$\mathfrak{A}\f$ is the effective domain of the conjugate of
    \f$\mathcal{R}\f$, i.e.,
    \f[
       \mathfrak{A} = \mathrm{dom}\,\mathcal{R}^*
                    = \{\vartheta\in\mathcal{X}^*\,:\,
                        \mathcal{R}^*(\vartheta) < \infty\}
    \f]
    where \f$\mathcal{R}^*\f$ denotes the Legendre-Fenchel transformation of
    \f$\mathcal{R}\f$.  This risk measure implements the penalized conditional
    value-at-risk
    \f[
       \mathcal{R}(X) = \sup_{\vartheta\in\mathfrak{A}}
         \left\{\mathbb{E}[\vartheta X]
          - \frac{\gamma}{2}\mathbb{E}[\vartheta^2]\right\}
    \f]
    for \f$\gamma > 0\f$.  This is implemented using the expectation risk
    quadrangle interface.  Thus, we represent \f$\mathcal{R}\f$ as
    \f[
       \mathcal{R}(X) = \inf_{t\in\mathbb{R}} \left\{
         t + \mathbb{E}\left[v(X-t)\right]
         \right\}
    \f]
    for an appropriately defined scalar regret function \f$v\f$.
    ROL implements this by augmenting the optimization vector \f$x_0\f$ with
    the parameter \f$t\f$, then minimizes jointly for \f$(x_0,t)\f$.
*/

namespace ROL {

template<class Real>
class MoreauYosidaCVaR : public ExpectationQuad<Real> {
private:

  Real prob_;
  Real eps_;

  Real omp_;
  Real ub_;

  void parseParameterList(ROL::ParameterList &parlist) {
    std::string type = parlist.sublist("SOL").get("Stochastic Component Type","Risk Averse");
    ROL::ParameterList list;
    if (type == "Risk Averse") {
      list = parlist.sublist("SOL").sublist("Risk Measure").sublist("Moreau-Yosida CVaR");
    }
    else if (type == "Error") {
      list = parlist.sublist("SOL").sublist("Error Measure").sublist("Moreau-Yosida-Koenker-Bassett");
    }
    else if (type == "Deviation") {
      list = parlist.sublist("SOL").sublist("Deviation Measure").sublist("Moreau-Yosida CVaR");
    }
    else if (type == "Regret") {
      list = parlist.sublist("SOL").sublist("Regret Measure").sublist("Moreau-Yosida Mean Absolute Loss");
    }
    prob_ = list.get<Real>("Confidence Level");
    eps_  = list.get<Real>("Smoothing Parameter");
  }

  void checkInputs(void) const {
    Real zero(0), one(1);
    ROL_TEST_FOR_EXCEPTION((prob_ <= zero) || (prob_ >= one), std::invalid_argument,
      ">>> ERROR (ROL::MoreauYosidaCVaR): Confidence level must be between 0 and 1!");
    ROL_TEST_FOR_EXCEPTION((eps_ <= zero), std::invalid_argument,
      ">>> ERROR (ROL::MoreauYosidaCVaR): Smoothing parameter must be positive!");
  }

  void setParameters(void) {
    Real one(1);
    omp_  = one-prob_;
    ub_   = eps_/omp_;
  }

public:
  /** \brief Constructor.

      @param[in]     prob    is the confidence level
      @param[in]     eps     is the regularization parameter
  */
  MoreauYosidaCVaR(Real prob, Real eps )
    : ExpectationQuad<Real>(), prob_(prob), eps_(eps) {
    checkInputs();
    setParameters();
  }

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"Moreau-Yosida CVaR" and
      within the "Moreau-Yosida CVaR" sublist should have the following parameters
      \li "Confidence Level" (between 0 and 1)
      \li "Smoothing Parameter" (must be positive)
  */
  MoreauYosidaCVaR(ROL::ParameterList &parlist)
    : ExpectationQuad<Real>() {
    parseParameterList(parlist);
    checkInputs();
    setParameters();
  }

  Real error(Real x, int deriv = 0) {
    Real zero(0), one(1);
    Real X = ((deriv==0) ? x : ((deriv==1) ? one : zero));
    return regret(x,deriv) - X;
  }

  Real regret(Real x, int deriv = 0) {
    Real zero(0), half(0.5), one(1), reg(0);
    if ( x <= zero ) {
      reg = 0;
    }
    else if ( x >= ub_ ) {
      reg  = ((deriv == 0) ? (x-half*ub_)/omp_
           : ((deriv == 1) ? one/omp_ : zero));
    }
    else {
      reg  = ((deriv == 0) ? half/eps_*x*x
           : ((deriv == 1) ? x/eps_ : one/eps_));
    }
    return reg;
  }

  void check(void) {
    ExpectationQuad<Real>::check();
    Real zero(0), one(1), two(2), p1(0.1);
    // Check v'(eps)
    Real x = eps_;
    Real vx = zero, vy = zero;
    Real dv = regret(x,1);
    Real t = one;
    Real diff = zero;
    Real err = zero;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v'(eps) is correct? \n";
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
    // check v''(eps) 
    vx = zero;
    vy = zero;
    dv = regret(x,2);
    t = one;
    diff = zero;
    err = zero;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v''(eps) is correct? \n";
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
    // check v''(eps) 
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
    // Check v'(0)
    x = -eps_;
    vx = zero;
    vy = zero;
    dv = regret(x,1);
    t = one;
    diff = zero;
    err = zero;
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
    // check v''(eps) 
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
  }

};

}
#endif
