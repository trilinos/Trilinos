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

#include "ROL_MixedQuantileQuadrangle.hpp"
#include "ROL_GaussLegendreQuadrature.hpp"

/** @ingroup risk_group
    \class ROL::SuperQuantileQuadrangle
    \brief Provides an interface for the risk measure associated with the
           super quantile quadrangle.

    The risk measure associated with the super quantile quadrangle is defined
    as
    \f[
       \mathcal{R}(X) = \frac{1}{1-\beta}\int_\beta^1\mathrm{CVaR}_{\alpha}(X)
          \,\mathrm{d}\alpha
    \f]
    where \f$0 \le \beta_n < 1\f$ and the conditional value-at-risk (CVaR) with
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
    Gauss-Legendre quadrature.  The corresponding quadrature points and weights
    are then used to construct a ROL::MixedQuantileQuadrangle risk measure.
    When using derivative-based optimization, the user can provide a smooth
    approximation of \f$(\cdot)_+\f$ using the ROL::PlusFunction class.
*/

namespace ROL {

template<class Real>
class SuperQuantileQuadrangle : public RiskMeasure<Real> {
private:
  Teuchos::RCP<MixedQuantileQuadrangle<Real> > mqq_;
  Teuchos::RCP<PlusFunction<Real> > plusFunction_;

  Real alpha_;
  int nQuad_;

  std::vector<Real> wts_;
  std::vector<Real> pts_;

  void checkInputs(void) const {
    TEUCHOS_TEST_FOR_EXCEPTION((alpha_ < 0 || alpha_ >= 1), std::invalid_argument,
      ">>> ERROR (ROL::SuperQuantileQuadrangle): Confidence level not between 0 and 1!");
    TEUCHOS_TEST_FOR_EXCEPTION(plusFunction_ == Teuchos::null, std::invalid_argument,
      ">>> ERROR (ROL::SuperQuantileQuadrangle): PlusFunction pointer is null!");
  }

  void initialize(void) {
     GaussLegendreQuadrature<Real> quad(nQuad_); // quad.test();
     quad.get(pts_,wts_);
     Real sum(0), half(0.5), one(1);
     for (int i = 0; i < nQuad_; ++i) {
       sum += wts_[i];
     }
     for (int i = 0; i < nQuad_; ++i) {
       wts_[i] /= sum;
       pts_[i] = one - alpha_*(half*(pts_[i] + one));
     }
     mqq_ = Teuchos::rcp(new MixedQuantileQuadrangle<Real>(pts_,wts_,plusFunction_));
  }

public:

  SuperQuantileQuadrangle( Teuchos::ParameterList &parlist )
    : RiskMeasure<Real>() {
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Super Quantile Quadrangle");
    // Grab confidence level and quadrature order
    alpha_ = list.get<Real>("Confidence Level");
    nQuad_ = list.get("Number of Quadrature Points",5);
    plusFunction_ = Teuchos::rcp(new PlusFunction<Real>(list));
    // Check inputs
    checkInputs();
    initialize();
  }

  SuperQuantileQuadrangle(const Real alpha,
                          const int nQuad,
                          const Teuchos::RCP<PlusFunction<Real> > &pf)
    : RiskMeasure<Real>(), plusFunction_(pf), alpha_(alpha), nQuad_(nQuad) {
    // Check inputs
    checkInputs();
    initialize();
  }

  Real computeStatistic(const Vector<Real> &x) const {
    std::vector<Real> xstat;
    Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic(xstat);
    Real stat(0);
    for (int i = 0; i < nQuad_; ++i) {
      stat += wts_[i] * xstat[i];
    }
    return stat;
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    mqq_->reset(x0,x);
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x,
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    mqq_->reset(x0,x,v0,v);
  }

  void update(const Real val, const Real weight) {
    mqq_->update(val,weight);
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    mqq_->update(val,g,weight);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
              const Real weight) {
    mqq_->update(val,g,gv,hv,weight);
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    return mqq_->getValue(sampler);
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    mqq_->getGradient(g,sampler);
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    mqq_->getHessVec(hv,sampler);
  }
};

}

#endif
