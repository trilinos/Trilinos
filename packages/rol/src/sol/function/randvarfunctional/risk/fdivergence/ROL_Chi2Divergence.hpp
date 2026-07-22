// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CHI2DIVERGENCE_HPP
#define ROL_CHI2DIVERGENCE_HPP

#include "ROL_FDivergence.hpp"

/** @ingroup risk_group
    \class ROL::Chi2Divergence
    \brief Provides an interface for the chi-squared-divergence distributionally robust
    expectation.

    This class defines a risk measure \f$\mathcal{R}\f$ that arises in distributionally
    robust stochastic programming.  \f$\mathcal{R}\f$ is given by
    \f[
       \mathcal{R}(X) = \sup_{\vartheta\in\mathfrak{A}}
           \mathbb{E}[\vartheta X]
    \f]
    where \f$\mathfrak{A}\f$ is called the ambiguity (or uncertainty) set and 
    is defined by a constraint on the \f$\chi^2\f$-divergence, i.e.,
    \f[
       \mathfrak{A} = \left\{\vartheta\in\mathcal{X}^*\,:\,
         \mathbb{E}[\vartheta] = 1,\; \vartheta \ge 0,\;\text{and}\;
         \frac{1}{2}\mathbb{E}[(\vartheta-1)^2] \le \epsilon\right\}.
    \f]
    \f$\mathcal{R}\f$ is a law-invariant, coherent risk measure.
*/

namespace ROL {

template<class Real>
class Chi2Divergence : public FDivergence<Real> {

public:
  /** \brief Constructor.

      @param[in]     thresh  is the tolerance for the F-divergence constraint
  */
  Chi2Divergence(const Real thresh) : FDivergence<Real>(thresh) {}

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"F-Divergence" and
      within the "F-Divergence" sublist should have the following parameters
      \li "Threshold" (greater than 0)
  */
  Chi2Divergence(ROL::ParameterList &parlist) : FDivergence<Real>(parlist) {}

  Real Fprimal(Real x, int deriv = 0) const {
    Real zero(0), one(1), half(0.5), val(0);
    if (deriv==0) {
      val = (x < zero) ? ROL_INF<Real>() : half*(x-one)*(x-one);
    }
    else if (deriv==1) {
      val = (x < zero) ? ROL_INF<Real>() : x-one;
    }
    else if (deriv==2) {
      val = (x < zero) ? ROL_INF<Real>() : one;
    }
    else {
      ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
        ">>> (ROL::Chi2Divergence): Derivative order must be 0, 1, or 2!");
    }
    return val;
  }

  Real Fdual(Real x, int deriv = 0) const {
    Real zero(0), one(1), half(0.5), val(0);
    if (deriv==0) {
      val = (x < -one) ? -half : (half*x + one)*x;
    }
    else if (deriv==1) {
      val = (x < -one) ? zero : x + one;
    }
    else if (deriv==2) {
      val = (x < -one) ? zero : one;
    }
    else {
      ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
        ">>> (ROL::Chi2Divergence): Derivative order must be 0, 1, or 2!");
    }
    return val;
  }
};

}

#endif
