// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_GENMOREAUYOSIDACVAR_HPP
#define ROL_GENMOREAUYOSIDACVAR_HPP

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
          - \frac{\gamma}{2}\mathbb{E}[(\vartheta-1)^2]\right\}
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
class GenMoreauYosidaCVaR : public ExpectationQuad<Real> {
private:

  Real prob_;
  Real lam_;
  Real eps_;

  Real alpha_;
  Real beta_;

  Real omp_;
  Real oma_;
  Real bmo_;
  Real lb_;
  Real ub_;

  void parseParameterList(ROL::ParameterList &parlist) {
    std::string type = parlist.sublist("SOL").get("Type","Risk Averse");
    ROL::ParameterList list;
    if (type == "Risk Averse") {
      list = parlist.sublist("SOL").sublist("Risk Measure").sublist("Generalized Moreau-Yosida CVaR");
    }
    else if (type == "Error") {
      list = parlist.sublist("SOL").sublist("Error Measure").sublist("Generalized Moreau-Yosida-Koenker-Bassett");
    }
    else if (type == "Deviation") {
      list = parlist.sublist("SOL").sublist("Deviation Measure").sublist("Generalized Moreau-Yosida CVaR");
    }
    else if (type == "Regret") {
      list = parlist.sublist("SOL").sublist("Regret Measure").sublist("Generalized Moreau-Yosida Mean Absolute Loss");
    }
    prob_ = list.get<Real>("Confidence Level");
    lam_  = list.get<Real>("Convex Combination Parameter");
    eps_  = list.get<Real>("Smoothing Parameter");
  }

  void checkInputs(void) const {
    Real zero(0), one(1);
    ROL_TEST_FOR_EXCEPTION((prob_ <= zero) || (prob_ >= one), std::invalid_argument,
      ">>> ERROR (ROL::GenMoreauYosidaCVaR): Confidence level must be between 0 and 1!");
    ROL_TEST_FOR_EXCEPTION((lam_ < zero) || (lam_ > one), std::invalid_argument,
      ">>> ERROR (ROL::GenMoreauYosidaCVaR): Convex combination parameter must be positive!");
    ROL_TEST_FOR_EXCEPTION((eps_ <= zero), std::invalid_argument,
      ">>> ERROR (ROL::GenMoreauYosidaCVaR): Smoothing parameter must be positive!");
  }

  void setParameters(void) {
    const Real one(1);
    omp_   = one-prob_;
    alpha_ = lam_;
    beta_  = (one-alpha_*prob_)/omp_;
    oma_   = one-alpha_;
    bmo_   = beta_-one;
    lb_    = -eps_*oma_;
    ub_    =  eps_*bmo_;
  }

public:
  /** \brief Constructor.

      @param[in]     prob    is the confidence level
      @param[in]     eps     is the regularization parameter
  */
  GenMoreauYosidaCVaR(Real prob, Real eps )
    : ExpectationQuad<Real>(), prob_(prob), lam_(0), eps_(eps) {
    checkInputs();
    setParameters();
  }

  /** \brief Constructor.

      @param[in]     prob    is the confidence level
      @param[in]     lam     is the convex combination parameter
      @param[in]     eps     is the regularization parameter
  */
  GenMoreauYosidaCVaR(Real prob, Real lam, Real eps )
    : ExpectationQuad<Real>(), prob_(prob), lam_(lam), eps_(eps) {
    checkInputs();
    setParameters();
  }

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"Moreau-Yosida CVaR" and
      within the "Moreau-Yosida CVaR" sublist should have the following parameters
      \li "Confidence Level" (between 0 and 1)
      \li "Convex Combination Parameter" (between 0 and 1)
      \li "Smoothing Parameter" (must be positive)
  */
  GenMoreauYosidaCVaR(ROL::ParameterList &parlist)
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
    if ( x <= lb_ ) {
      reg = ((deriv == 0) ? alpha_*x + half*lb_*oma_
          : ((deriv == 1) ? alpha_ : zero));
    }
    else if ( x >= ub_ ) {
      reg = ((deriv == 0) ? beta_*x - half*ub_*bmo_
          : ((deriv == 1) ? beta_ : zero));
    }
    else {
      reg = ((deriv == 0) ? half/eps_*x*x + x
          : ((deriv == 1) ? x/eps_ + one : one/eps_));
    }
    return reg;
  }

  void check(void) {
    ExpectationQuad<Real>::check();
    Real zero(0), one(1), two(2), p1(0.1);
    // Check v'(ub)
    Real x = ub_;
    Real vx = zero, vy = zero;
    Real dv = regret(x,1);
    Real t = one;
    Real diff = zero;
    Real err = zero;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v'(ub) is correct? \n";
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
    // check v''(ub) 
    vx = zero;
    vy = zero;
    dv = regret(x,2);
    t = one;
    diff = zero;
    err = zero;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v''(ub) is correct? \n";
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
    // Check v'(lb)
    x = lb_;
    vx = zero;
    vy = zero;
    dv = regret(x,1);
    t = one;
    diff = zero;
    err = zero;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v'(lb) is correct? \n";
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
    // check v''(lb) 
    vx = zero;
    vy = zero;
    dv = regret(x,2);
    t = one;
    diff = zero;
    err = zero;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v''(lb) is correct? \n";
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
