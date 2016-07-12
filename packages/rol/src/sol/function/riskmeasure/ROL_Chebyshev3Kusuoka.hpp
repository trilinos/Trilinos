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

#ifndef ROL_CHEBYSHEV3KUSUOKA_HPP
#define ROL_CHEBYSHEV3KUSUOKA_HPP

#include "ROL_SingletonKusuoka.hpp"
#include "ROL_GaussChebyshev3Quadrature.hpp"

/** @ingroup risk_group
    \class ROL::Chebyshev3Kusuoka
    \brief Provides an interface for the Chebyshev 3 Kusuoka risk measure.

    The Chebyshev 3 Kusuoka risk measure is defined as
    \f[
       \mathcal{R}(X) = \int_0^1 w(\alpha) \mathrm{CVaR}_{\alpha}(X)
          \,\mathrm{d}\alpha
    \f]
    where the conditional value-at-risk (CVaR) with confidence level
    \f$0\le \alpha < 1\f$ is
    \f[
       \mathrm{CVaR}_\alpha(X) = \inf_{t\in\mathbb{R}} \left\{
         t + \frac{1}{1-\alpha} \mathbb{E}\left[(X-t)_+\right]
         \right\}, \quad (x)_+ = \max\{0,x\},
    \f]
    and the weight function \f$w\f$ is
    \f[
       w(x) = \sqrt{\frac{x}{1-x}}.
    \f]
    If the distribution of \f$X\f$ is continuous, then
    \f$\mathrm{CVaR}_{\alpha}(X)\f$ is the conditional
    expectation of \f$X\f$ exceeding the \f$\alpha\f$-quantile of \f$X\f$ and
    the optimal \f$t\f$ is the \f$\alpha\f$-quantile.
    Additionally, \f$\mathcal{R}\f$ is a law-invariant coherent risk measure.

    ROL implements \f$\mathcal{R}\f$ by approximating the integral with
    Gauss-Chebyshev quadrature of the third kind.  The corresponding
    quadrature points and weights are then used to construct a
    ROL::MixedQuantileQuadrangle risk measure.
    When using derivative-based optimization, the user can provide a smooth
    approximation of \f$(\cdot)_+\f$ using the ROL::PlusFunction class.
*/

namespace ROL {

template<class Real>
class Chebyshev3Kusuoka : public SingletonKusuoka<Real> {
private:
  Teuchos::RCP<PlusFunction<Real> > plusFunction_;

  Real lower_, upper_;
  int nQuad_; 

  std::vector<Real> wts_;
  std::vector<Real> pts_;

  void checkInputs(void) const {
    TEUCHOS_TEST_FOR_EXCEPTION(lower_ > upper_, std::invalid_argument,
      ">>> ERROR (ROL::Chebyshev3Kusuoka): Lower bound exceeds upper!");
    TEUCHOS_TEST_FOR_EXCEPTION(lower_ < static_cast<Real>(0), std::invalid_argument,
      ">>> ERROR (ROL::Chebyshev3Kusuoka): Lower bound is less than zero!");
    TEUCHOS_TEST_FOR_EXCEPTION(static_cast<Real>(1) < upper_, std::invalid_argument,
      ">>> ERROR (ROL::Chebyshev3Kusuoka): Upper bound is greater than one!");
    TEUCHOS_TEST_FOR_EXCEPTION(plusFunction_ == Teuchos::null, std::invalid_argument,
      ">>> ERROR (ROL::Chebyshev3Kusuoka): PlusFunction pointer is null!");
  }

  void initialize(void) {
     GaussChebyshev3Quadrature<Real> quad(nQuad_); // quad.test();
     quad.get(pts_,wts_);
     Real sum(0), half(0.5), one(1);
     for (int i = 0; i < nQuad_; ++i) {
       sum += wts_[i];
     }
     for (int i = 0; i < nQuad_; ++i) {
       wts_[i] /= sum;
       pts_[i] = lower_ + (upper_-lower_)*half*(pts_[i] + one);
     }
     SingletonKusuoka<Real>::buildMixedQuantile(pts_,wts_,plusFunction_);
  }

public:
  Chebyshev3Kusuoka( Teuchos::ParameterList &parlist )
    : SingletonKusuoka<Real>() {
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Chebyshev 3 Kusuoka");
    // Grab confidence level and quadrature order
    lower_ = list.get("Lower Bound",0.0);
    upper_ = list.get("Upper Bound",1.0);
    nQuad_ = list.get("Number of Quadrature Points",5);
    plusFunction_ = Teuchos::rcp(new PlusFunction<Real>(list));
    // Check inputs
    checkInputs();
    initialize();
  }

  Chebyshev3Kusuoka(const Real lower, const Real upper,
                    const int nQuad,
                    const Teuchos::RCP<PlusFunction<Real> > &pf)
    : RiskMeasure<Real>(), plusFunction_(pf), lower_(lower), upper_(upper), nQuad_(nQuad) {
    // Check inputs
    checkInputs();
    initialize();
  }
};

}

#endif
