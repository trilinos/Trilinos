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

#ifndef ROL_SINGLETONKUSUOKA_HPP
#define ROL_SINGLETONKUSUOKA_HPP

#include "ROL_MixedQuantileQuadrangle.hpp"

/** @ingroup risk_group
    \class ROL::SingletonKusuoka
    \brief Provides an interface for singleton Kusuoka risk measures.

    Kusuoka's representation for law-invariant risk measures is
    \f[
       \mathcal{R}(X) = \sup_{\mu\in\mathfrak{M}}
          \int_0^1 \mathrm{CVaR}_{\alpha}(X)\,\mathrm{d}\mu(\alpha)
    \f]
    where the conditional value-at-risk (CVaR) with confidence level
    \f$0\le \alpha < 1\f$ is
    \f[
       \mathrm{CVaR}_\alpha(X) = \inf_{t\in\mathbb{R}} \left\{
         t + \frac{1}{1-\alpha} \mathbb{E}\left[(X-t)_+\right]
         \right\}, \quad (x)_+ = \max\{0,x\},
    \f]
    and \f$\mathfrak{M}\f$ is a subset of distributions on the interval
    \f$[0,1)\f$.  By singleton Kusuoka, we refer to the case where the set
    \f$\mathfrak{M}\f$ is a singleton.  If the distribution
    \f$\mu\in\mathfrak{M}\f$ is discrete, then the corresponding risk measure
    is a mixed quantile quadrangle risk measure.

    If the distribution of \f$X\f$ is continuous, then
    \f$\mathrm{CVaR}_{\alpha}(X)\f$ is the conditional
    expectation of \f$X\f$ exceeding the \f$\alpha\f$-quantile of \f$X\f$ and
    the optimal \f$t\f$ is the \f$\alpha\f$-quantile.
    Additionally, \f$\mathcal{R}\f$ is a law-invariant coherent risk measure.

    ROL implements \f$\mathcal{R}\f$ by approximating the integral with
    Gauss-Chebyshev quadrature of the first kind.  The corresponding quadrature
    points and weights are then used to construct a
    ROL::MixedQuantileQuadrangle risk measure.
    When using derivative-based optimization, the user can provide a smooth
    approximation of \f$(\cdot)_+\f$ using the ROL::PlusFunction class.
*/

namespace ROL {

template<class Real>
class SingletonKusuoka : public RiskMeasure<Real> {
private:
  Teuchos::RCP<MixedQuantileQuadrangle<Real> > mqq_;
  Teuchos::RCP<PlusFunction<Real> > plusFunction_;

  std::vector<Real> wts_;
  std::vector<Real> pts_;

  void checkInputs(void) const {
    TEUCHOS_TEST_FOR_EXCEPTION(plusFunction_ == Teuchos::null, std::invalid_argument,
      ">>> ERROR (ROL::SingletonKusuoka): PlusFunction pointer is null!");
  }

protected:
  void buildMixedQuantile(const std::vector<Real> &pts, const std::vector<Real> &wts,
                          const Teuchos::RCP<PlusFunction<Real> > &pf) {
     pts_.clear(); pts_.assign(pts.begin(),pts.end());
     wts_.clear(); wts_.assign(wts.begin(),wts.end());
     plusFunction_ = pf;
     mqq_ = Teuchos::rcp(new MixedQuantileQuadrangle<Real>(pts,wts,pf));
  }

public:
  SingletonKusuoka(void) : RiskMeasure<Real>() {}

  SingletonKusuoka( const std::vector<Real> &pts, const std::vector<Real> &wts,
                    const Teuchos::RCP<PlusFunction<Real> > &pf)
    : RiskMeasure<Real>() {
    // Check inputs
    checkInputs();
    buildMixedQuantile(pts,wts,pf);
  }

  Real computeStatistic(const Vector<Real> &x) const {
    std::vector<Real> xstat;
    Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic(xstat);
    Real stat(0);
    int nQuad = static_cast<int>(wts_.size());
    for (int i = 0; i < nQuad; ++i) {
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
