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

#ifndef ROL_SUPERQUANTILEQUADRANGLE_HPP
#define ROL_SUPERQUANTILEQUADRANGLE_HPP

#include "ROL_SpectralRisk.hpp"
#include "ROL_GaussLegendreQuadrature.hpp"
#include "ROL_Fejer2Quadrature.hpp"

/** @ingroup risk_group
    \class ROL::SecondOrderCVaR
    \brief Provides an interface for the risk measure associated with the
           super quantile quadrangle.

    The risk measure associated with the super quantile quadrangle is defined
    as
    \f[
       \mathcal{R}(X) = \frac{1}{1-\beta}\int_\beta^1\mathrm{CVaR}_{\alpha}(X)
          \,\mathrm{d}\alpha
    \f]
    where \f$0 \le \beta < 1\f$ and the conditional value-at-risk (CVaR) with
    confidence level \f$0\le \alpha < 1\f$ is
    \f[
       \mathrm{CVaR}_\alpha(X) = \inf_{t\in\mathbb{R}} \left\{
         t + \frac{1}{1-\alpha} \mathbb{E}\left[(X-t)_+\right]
         \right\}
    \f]
    where \f$(x)_+ = \max\{0,x\}\f$.  If the distribution of \f$X\f$ is
    continuous, then \f$\mathrm{CVaR}_{\alpha}(X)\f$ is the conditional
    expectation of \f$X\f$ exceeding the \f$\alpha\f$-quantile of \f$X\f$ and
    the optimal \f$t\f$ is the \f$\alpha\f$-quantile.
    Additionally, \f$\mathcal{R}\f$ is a law-invariant coherent risk measure.

    ROL implements \f$\mathcal{R}\f$ by approximating the integral with
    Gauss-Legendre or Fejer 2 quadrature.  The corresponding quadrature points
    and weights are then used to construct a ROL::MixedQuantileQuadrangle risk
    measure.  When using derivative-based optimization, the user can provide a
    smooth approximation of \f$(\cdot)_+\f$ using the ROL::PlusFunction class.
*/

namespace ROL {

template<class Real>
class SecondOrderCVaR : public SpectralRisk<Real> {
private:
  ROL::Ptr<PlusFunction<Real> > plusFunction_;

  Real alpha_;
  int nQuad_;
  bool useGauss_;

  std::vector<Real> wts_;
  std::vector<Real> pts_;

  void checkInputs(void) const {
    ROL_TEST_FOR_EXCEPTION((alpha_ < 0 || alpha_ >= 1), std::invalid_argument,
      ">>> ERROR (ROL::SecondOrderCVaR): Confidence level not between 0 and 1!");
    ROL_TEST_FOR_EXCEPTION(plusFunction_ == ROL::nullPtr, std::invalid_argument,
      ">>> ERROR (ROL::SecondOrderCVaR): PlusFunction pointer is null!");
  }

  void initializeQuad(void) {
    ROL::Ptr<Quadrature1D<Real> > quad;
    if ( useGauss_ ) {
      quad = ROL::makePtr<GaussLegendreQuadrature<Real>>(nQuad_);
    }
    else {
      quad = ROL::makePtr<Fejer2Quadrature<Real>>(nQuad_);
    }
    // quad->test();
    quad->get(pts_,wts_);
    Real sum(0), half(0.5), one(1);
    for (int i = 0; i < nQuad_; ++i) {
      sum += wts_[i];
    }
    for (int i = 0; i < nQuad_; ++i) {
      wts_[i] /= sum;
      pts_[i] = one - alpha_*(half*(pts_[i] + one));
    }
    SpectralRisk<Real>::buildMixedQuantile(pts_,wts_,plusFunction_);
  }

public:

  SecondOrderCVaR( ROL::ParameterList &parlist )
    : SpectralRisk<Real>() {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Second Order CVaR");
    // Grab confidence level and quadrature order
    alpha_ = list.get<Real>("Confidence Level");
    nQuad_ = list.get("Number of Quadrature Points",5);
    useGauss_ = list.get("Use Gauss-Legendre Quadrature",true);
    plusFunction_ = ROL::makePtr<PlusFunction<Real>>(list);
    // Check inputs
    checkInputs();
    initializeQuad();
  }

  SecondOrderCVaR(const Real alpha,
                          const int nQuad,
                          const ROL::Ptr<PlusFunction<Real> > &pf,
                          const bool useGauss = true)
    : SpectralRisk<Real>(), plusFunction_(pf),
      alpha_(alpha), nQuad_(nQuad), useGauss_(useGauss) {
    // Check inputs
    checkInputs();
    initializeQuad();
  }
};

}

#endif
