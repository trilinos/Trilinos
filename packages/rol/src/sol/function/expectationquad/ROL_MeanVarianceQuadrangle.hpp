// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MEANVARIANCEQUAD_HPP
#define ROL_MEANVARIANCEQUAD_HPP

#include "ROL_ExpectationQuad.hpp"

/** @ingroup risk_group
    \class ROL::MeanVarianceQuadrangle
    \brief Provides an interface for the mean plus variance risk measure
           using the expectation risk quadrangle.

    The mean plus variances risk measure is
    \f[
       \mathcal{R}(X) = \mathbb{E}[X]
        + c \mathbb{E}[|X-\mathbb{E}[X]|^2]
    \f]
    where \f$c \ge 0\f$.
    \f$\mathcal{R}\f$ is law-invariant, but not coherent since it
    violates positive homogeneity.  The associated scalar regret
    function is
    \f[
       v(x) = c x^2 + x
    \f]
    and the mean-plus-variance risk measure is computed as
    \f[
       \mathcal{R}(X) = \inf_{t\in\mathbb{R}}\left\{
           t + \mathbb{E}[v(X-t)] \right\}.
    \f]
    ROL implements this by augmenting the optimization vector \f$x_0\f$ with
    the parameter \f$t\f$, then minimizes jointly for \f$(x_0,t)\f$.
*/

namespace ROL {

template<class Real>
class MeanVarianceQuadrangle : public ExpectationQuad<Real> {
private:
  Real coeff_;

  void parseParameterList(ROL::ParameterList &parlist) {
    std::string type = parlist.sublist("SOL").get("Type","Risk Averse");
    ROL::ParameterList list;
    if (type == "Risk Averse") {
      list = parlist.sublist("SOL").sublist("Risk Measure").sublist("Safety Margin");
    }
    else if (type == "Regret") {
      list = parlist.sublist("SOL").sublist("Regret Measure").sublist("Mean L2");
    }
    else if (type == "Error" || type == "Deviation") {
      coeff_ = static_cast<Real>(1);
      return;
    }
    coeff_ = list.get<Real>("Coefficient");
  }

  void checkInputs(void) const {
    Real zero(0);
    ROL_TEST_FOR_EXCEPTION((coeff_ <= zero), std::invalid_argument,
      ">>> ERROR (ROL::MeanVarianceQuadrangle): Coefficient must be positive!");
  }

public:
  /** \brief Constructor.

      @param[in]     coeff   is the weight for variance term
  */
  MeanVarianceQuadrangle(const Real coeff = 1)
    : ExpectationQuad<Real>(), coeff_(coeff) {
    checkInputs();
  }

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"Mean-Variance Quadrangle" and
      within the "Mean-Variance Quadrangle" sublist should have the following parameters
      \li "Coefficient" (array of positive scalars).
  */
  MeanVarianceQuadrangle(ROL::ParameterList &parlist)
    : ExpectationQuad<Real>() {
    parseParameterList(parlist);
    checkInputs();
  }

  Real error(Real x, int deriv = 0) {
    Real err(0), two(2);
    if (deriv==0) {
      err = coeff_*x*x;
    }
    else if (deriv==1) {
      err = two*coeff_*x;
    }
    else {
      err = two*coeff_;
    }
    return err;
  }

  Real regret(Real x, int deriv = 0) {
    Real zero(0), one(1);
    Real X = ((deriv==0) ? x : ((deriv==1) ? one : zero));
    Real reg = error(x,deriv) + X;
    return reg;
  }

};

}
#endif
