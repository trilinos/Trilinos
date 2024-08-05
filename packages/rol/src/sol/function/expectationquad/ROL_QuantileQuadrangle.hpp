// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SMOOTHCVARQUAD_HPP
#define ROL_SMOOTHCVARQUAD_HPP

#include "ROL_ExpectationQuad.hpp"
#include "ROL_PlusFunction.hpp"

/** @ingroup risk_group
    \class ROL::QuantileQuadrangle
    \brief Provides an interface for a convex combination of the
           expected value and the conditional value-at-risk using the
           expectation risk quadrangle.

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

    This class defines a convex combination of expected value and the conditional
    value-at-risk using the expectation risk quadrangle.  In this case, the scalar
    regret function is
    \f[
       v(x) = \alpha (x)_+ - \lambda (-x)_+
    \f]
    for \f$\alpha > 1\f$ and \f$0 \le \lambda < 1\f$.  The associated confidence
    level for the conditional value-at-risk is
    \f[
       \beta = \frac{\alpha-1}{\alpha-\lambda}.
    \f]
    This convex combination of expected value and the conditional value-at-risk
    is then realized as
    \f[
       \mathcal{R}(X) = \inf_{t\in\mathbb{R}}\left\{
           t + \mathbb{E}[v(X-t)] \right\}.
    \f]
    ROL implements this by augmenting the optimization vector \f$x_0\f$ with
    the parameter \f$t\f$, then minimizes jointly for \f$(x_0,t)\f$.

    When using derivative-based optimization, the user can provide a smooth
    approximation of \f$(\cdot)_+\f$ using the ROL::PlusFunction class.
*/

namespace ROL {

template<class Real>
class QuantileQuadrangle : public ExpectationQuad<Real> {
private:

  ROL::Ptr<PlusFunction<Real> > pf_;

  Real prob_;
  Real lam_;
  Real eps_;

  Real alpha_;
  Real beta_;

  void parseParameterList(ROL::ParameterList &parlist) {
    std::string type = parlist.sublist("SOL").get("Type","Risk Averse");
    ROL::ParameterList list;
    if (type == "Risk Averse") {
      list = parlist.sublist("SOL").sublist("Risk Measure").sublist("CVaR");
    }
    else if (type == "Error") {
      list = parlist.sublist("SOL").sublist("Error Measure").sublist("Koenker-Bassett");
    }
    else if (type == "Deviation") {
      list = parlist.sublist("SOL").sublist("Deviation Measure").sublist("CVaR");
    }
    else if (type == "Regret") {
      list = parlist.sublist("SOL").sublist("Regret Measure").sublist("Mean Absolute Loss");
    }
    // Get CVaR parameters
    prob_ = list.get<Real>("Confidence Level");
    lam_  = list.get<Real>("Convex Combination Parameter");
    eps_  = list.get<Real>("Smoothing Parameter");
    // Build plus function
    pf_   = ROL::makePtr<PlusFunction<Real>>(list);
  }

  void checkInputs(void) const {
    Real zero(0), one(1);
    ROL_TEST_FOR_EXCEPTION((prob_ <= zero) || (prob_ >= one), std::invalid_argument,
      ">>> ERROR (ROL::QuantileQuadrangle): Confidence level must be between 0 and 1!");
    ROL_TEST_FOR_EXCEPTION((lam_ < zero) || (lam_ > one), std::invalid_argument,
      ">>> ERROR (ROL::QuantileQuadrangle): Convex combination parameter must be positive!");
    ROL_TEST_FOR_EXCEPTION((eps_ <= zero), std::invalid_argument,
      ">>> ERROR (ROL::QuantileQuadrangle): Smoothing parameter must be positive!");
    ROL_TEST_FOR_EXCEPTION(pf_ == ROL::nullPtr, std::invalid_argument,
      ">>> ERROR (ROL::QuantileQuadrangle): PlusFunction pointer is null!");
  }

  void setParameters(void) {
    Real one(1);
    alpha_ = lam_;
    beta_  = (one-alpha_*prob_)/(one-prob_);
  }

public:
  /** \brief Constructor.

      @param[in]     prob    is the confidence level
      @param[in]     eps     is the smoothing parameter for the plus function approximation
      @param[in]     pf      is the plus function or an approximation
  */
  QuantileQuadrangle(Real prob, Real eps, ROL::Ptr<PlusFunction<Real> > &pf ) 
    : ExpectationQuad<Real>(), prob_(prob), lam_(0), eps_(eps), pf_(pf) {
    checkInputs();
    setParameters();
  }

  /** \brief Constructor.

      @param[in]     prob    is the confidence level
      @param[in]     lam     is the convex combination parameter (coeff=0
                             corresponds to the expected value whereas coeff=1
                             corresponds to the conditional value-at-risk)
      @param[in]     eps     is the smoothing parameter for the plus function approximation
      @param[in]     pf      is the plus function or an approximation
  */
  QuantileQuadrangle(Real prob, Real lam, Real eps,
                     ROL::Ptr<PlusFunction<Real> > &pf ) 
    : ExpectationQuad<Real>(), prob_(prob), lam_(lam), eps_(eps), pf_(pf) {
    checkInputs();
    setParameters();
  }

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"CVaR" and
      within the "CVaR" sublist should have the following parameters
      \li "Confidence Level" (between 0 and 1)
      \li "Convex Combination Parameter" (between 0 and 1)
      \li "Smoothing Parameter" (must be positive)
      \li A sublist for plus function information.
  */
  QuantileQuadrangle(ROL::ParameterList &parlist) : ExpectationQuad<Real>() {
    parseParameterList(parlist);
    // Check inputs
    checkInputs();
    setParameters();
  }

  Real error(Real x, int deriv = 0) {
    Real one(1);
    Real err = (beta_-one)*pf_->evaluate(x,deriv) 
              + ((deriv%2) ? -one : one)*(one-alpha_)*pf_->evaluate(-x,deriv);
    return err;
  }

  Real regret(Real x, int deriv = 0) {
    Real zero(0), one(1);
    Real X = ((deriv==0) ? x : ((deriv==1) ? one : zero));
    Real reg = error(x,deriv) + X;
    return reg;
  }

  void check(void) {
    ExpectationQuad<Real>::check();
    // Check v'(eps)
    Real x = eps_, two(2), p1(0.1), one(1), zero(0);
    Real vx(0), vy(0);
    Real dv = regret(x,1);
    Real t(1), diff(0), err(0);
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
