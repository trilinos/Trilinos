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

  Real Fprimal(Real x, int deriv = 0) {
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

  Real Fdual(Real x, int deriv = 0) {
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
