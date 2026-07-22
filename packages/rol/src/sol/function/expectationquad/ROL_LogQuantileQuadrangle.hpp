// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LOGQUANTILEQUAD_HPP
#define ROL_LOGQUANTILEQUAD_HPP

#include "ROL_ExpectationQuad.hpp"
#include "ROL_PlusFunction.hpp"

/** @ingroup risk_group
    \class ROL::LogQuantileQuadrangle
    \brief Provides an interface for the conditioanl entropic risk using
           the expectation risk quadrangle.

    This class defines the conditional entropic risk measure using the
    framework of the expectation risk quadrangle.  In this case, the scalar
    regret function is
    \f[
       v(x) = \lambda(\exp\left(\frac{x}{\lambda}\right)-1)_+ - \alpha (-x)_+
    \f]
    for \f$\lambda > 0\f$ and \f$0 \le \alpha < 1\f$. 
    The entropic risk measure is then implemented as
    \f[
       \mathcal{R}(X) = \inf_{t\in\mathbb{R}}\left\{
           t + \mathbb{E}[v(X-t)] \right\}.
    \f]
    The conditional entropic risk is convex, translation equivariant and
    monotonic.
    ROL implements this by augmenting the optimization vector \f$x_0\f$ with
    the parameter \f$t\f$, then minimizes jointly for \f$(x_0,t)\f$.
*/

namespace ROL {

template<class Real>
class LogQuantileQuadrangle : public ExpectationQuad<Real> {
private:

  ROL::Ptr<PlusFunction<Real> > pf_;

  Real alpha_;
  Real rate_;
  Real eps_;

  void parseParameterList(ROL::ParameterList &parlist) {
    std::string type = parlist.sublist("SOL").get("Type","Risk Averse");
    ROL::ParameterList list;
    if (type == "Risk Averse") {
      list = parlist.sublist("SOL").sublist("Risk Measure").sublist("Log Quantile");
    }
    else if (type == "Error") {
      list = parlist.sublist("SOL").sublist("Error Measure").sublist("Log Quantile");
    }
    else if (type == "Deviation") {
      list = parlist.sublist("SOL").sublist("Deviation Measure").sublist("Log Quantile");
    }
    else if (type == "Regret") {
      list = parlist.sublist("SOL").sublist("Regret Measure").sublist("Log Quantile");
    }
    // Check CVaR inputs
    alpha_  = list.get<Real>("Slope for Linear Growth");
    rate_   = list.get<Real>("Rate for Exponential Growth");
    eps_    = list.get<Real>("Smoothing Parameter");
    // Build plus function
    pf_ = ROL::makePtr<PlusFunction<Real>>(list);
  }

  void checkInputs(void) const {
    Real zero(0), one(1);
    ROL_TEST_FOR_EXCEPTION((alpha_ < zero) || (alpha_ >= one), std::invalid_argument,
      ">>> ERROR (ROL::LogQuantileQuadrangle): Linear growth rate must be between 0 and 1!");
    ROL_TEST_FOR_EXCEPTION((rate_ <= zero), std::invalid_argument,
      ">>> ERROR (ROL::LogQuantileQuadrangle): Exponential growth rate must be positive!");
    ROL_TEST_FOR_EXCEPTION((eps_ <= zero), std::invalid_argument,
      ">>> ERROR (ROL::LogQuantileQuadrangle): Smoothing parameter must be positive!");
    ROL_TEST_FOR_EXCEPTION(pf_ == ROL::nullPtr, std::invalid_argument,
      ">>> ERROR (ROL::LogQuantileQuadrangle): PlusFunction pointer is null!");
  }

public:
  /** \brief Constructor.

      @param[in]     alpha    is the scale parameter for the negative branch of the scalar regret
      @param[in]     rate     is the rate parameter for the positive branch of the scalar regret
      @param[in]     eps      is the smoothing parameter for the approximate plus function
      @param[in]     pf       is the plus function or an approximation
  */
  LogQuantileQuadrangle(Real alpha, Real rate, Real eps,
                        ROL::Ptr<PlusFunction<Real> > &pf ) 
    : ExpectationQuad<Real>(), alpha_(alpha), rate_(rate), eps_(eps), pf_(pf) {
    checkInputs();
  }

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measures"->"Log-Quantile Quadrangle"
      and withing the "Log-Quantile Quadrangle" sublist should have
      \li "Slope for Linear Growth" (between 0 and 1)
      \li "Rate for Exponential Growth" (must be positive)
      \li "Smoothing Parameter" (must be positive)
      \li A sublist for plus function information.
  */
  LogQuantileQuadrangle(ROL::ParameterList &parlist)
    : ExpectationQuad<Real>() {
    parseParameterList(parlist);
    checkInputs();
  }

  Real error(Real x, int deriv = 0) {
    Real zero(0), one(1);
    ROL_TEST_FOR_EXCEPTION( (deriv > 2), std::invalid_argument,
      ">>> ERROR (ROL::LogQuantileQuadrangle::error): deriv greater than 2!");
    ROL_TEST_FOR_EXCEPTION( (deriv < 0), std::invalid_argument,
      ">>> ERROR (ROL::LogQuantileQuadrangle::error): deriv less than 0!");

    Real X = ((deriv == 0) ? x : ((deriv == 1) ? one : zero));
    return regret(x,deriv) - X;
  }

  Real regret(Real x, int deriv = 0) {
    Real zero(0), one(1);
    ROL_TEST_FOR_EXCEPTION( (deriv > 2), std::invalid_argument,
      ">>> ERROR (ROL::LogQuantileQuadrangle::regret): deriv greater than 2!");
    ROL_TEST_FOR_EXCEPTION( (deriv < 0), std::invalid_argument,
      ">>> ERROR (ROL::LogQuantileQuadrangle::regret): deriv less than 0!");

    Real arg  = std::exp(rate_*x);
    Real sarg = rate_*arg;
    Real reg  = (pf_->evaluate(arg-one,deriv) *
                  ((deriv == 0) ? one/rate_ : ((deriv == 1) ? arg : sarg*arg))
                + ((deriv == 2) ? pf_->evaluate(arg-one,deriv-1)*sarg : zero))
                + ((deriv%2 == 0) ? -one : one) * alpha_ * pf_->evaluate(-x,deriv);
    return reg;
  }

  void check(void) {
    ExpectationQuad<Real>::check();
    // Check v'(eps)
    Real x = eps_, two(2), p1(0.1), zero(0), one(1);
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
