// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TRUNCATEDMEANQUAD_HPP
#define ROL_TRUNCATEDMEANQUAD_HPP

#include "ROL_ExpectationQuad.hpp"

namespace ROL {

template<class Real>
class TruncatedMeanQuadrangle : public ExpectationQuad<Real> {
private:

  Real beta_;

  void parseParameterList(ROL::ParameterList &parlist) {
    std::string type = parlist.sublist("SOL").get("Type","Risk Averse");
    ROL::ParameterList list;
    if (type == "Risk Averse") {
      list = parlist.sublist("SOL").sublist("Risk Measure").sublist("Truncated Mean");
    }
    else if (type == "Error") {
      list = parlist.sublist("SOL").sublist("Error Measure").sublist("Huber");
    }
    else if (type == "Deviation") {
      list = parlist.sublist("SOL").sublist("Deviation Measure").sublist("Truncated Mean");
    }
    else if (type == "Regret") {
      list = parlist.sublist("SOL").sublist("Regret Measure").sublist("Truncated Mean");
    }
    beta_ = list.get<Real>("Threshold");
  }

  void checkInputs(void) const {
    Real zero(0);
    ROL_TEST_FOR_EXCEPTION((beta_ <= zero), std::invalid_argument,
      ">>> ERROR (ROL::TruncatedMeanQuadrangle): Threshold must be positive!");
  }

public:

  TruncatedMeanQuadrangle(Real beta)
    : ExpectationQuad<Real>(), beta_(beta) {
    checkInputs();
  }

  TruncatedMeanQuadrangle(ROL::ParameterList &parlist)
    : ExpectationQuad<Real>() {
    parseParameterList(parlist);
    checkInputs();
  }

  Real error(Real x, int deriv = 0) {
    bool inside = ( std::abs(x) < beta_ ? true : false );
    Real err(0), zero(0), half(0.5), one(1), two(2);
    if (deriv==0) {
      err = (inside ? half*std::pow(x,two)/beta_ : std::abs(x)-half*beta_);
    }
    else if (deriv==1) {
      err = (inside ? x/beta_ : ((zero < x) - (x < zero)));
    }
    else {
      err = (inside ? one/beta_ : zero);
    }
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
    // Check v'(beta)
    Real x = beta_, zero(0), one(1), two(2), p1(0.1);
    Real vx = zero, vy = zero;
    Real dv = regret(x,1);
    Real t = one;
    Real diff = zero;
    Real err = zero;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v'(beta) is correct? \n";
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
    // Check v'(-beta)
    x = -beta_;
    vx = zero;
    vy = zero;
    dv = regret(x,1);
    t = one;
    diff = zero;
    err = zero;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v'(-beta) is correct? \n";
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
  }

};

}
#endif
