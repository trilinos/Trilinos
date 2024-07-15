// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SPECTRALRISK_HPP
#define ROL_SPECTRALRISK_HPP

#include "ROL_MixedCVaR.hpp"
#include "ROL_DistributionFactory.hpp"

/** @ingroup risk_group
    \class ROL::SpectralRisk
    \brief Provides an interface for spectral risk measures.

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
    \f$[0,1)\f$.  By spectral risk measures, we refer to the case where the set
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
    ROL::MixedCVaR risk measure.
    When using derivative-based optimization, the user can provide a smooth
    approximation of \f$(\cdot)_+\f$ using the ROL::PlusFunction class.
*/

namespace ROL {

template<class Real>
class SpectralRisk : public RandVarFunctional<Real> {
private:
  ROL::Ptr<MixedCVaR<Real> > mqq_;
  ROL::Ptr<PlusFunction<Real> > plusFunction_;

  std::vector<Real> wts_;
  std::vector<Real> pts_;

  void checkInputs(const ROL::Ptr<Distribution<Real> > &dist = ROL::nullPtr) const {
    ROL_TEST_FOR_EXCEPTION(plusFunction_ == ROL::nullPtr, std::invalid_argument,
      ">>> ERROR (ROL::SpectralRisk): PlusFunction pointer is null!");
    if ( dist != ROL::nullPtr) {
      Real lb = dist->lowerBound();
      Real ub = dist->upperBound();
      ROL_TEST_FOR_EXCEPTION(lb < static_cast<Real>(0), std::invalid_argument,
        ">>> ERROR (ROL::SpectralRisk): Distribution lower bound less than zero!");
      ROL_TEST_FOR_EXCEPTION(ub > static_cast<Real>(1), std::invalid_argument,
        ">>> ERROR (ROL::SpectralRisk): Distribution upper bound greater than one!");
    }
  }

protected:
  void buildMixedQuantile(const std::vector<Real> &pts, const std::vector<Real> &wts,
                          const ROL::Ptr<PlusFunction<Real> > &pf) {
     pts_.clear(); pts_.assign(pts.begin(),pts.end());
     wts_.clear(); wts_.assign(wts.begin(),wts.end());
     plusFunction_ = pf;
     mqq_ = ROL::makePtr<MixedCVaR<Real>>(pts,wts,pf);
  }

  void buildQuadFromDist(std::vector<Real> &pts, std::vector<Real> &wts,
                   const int nQuad, const ROL::Ptr<Distribution<Real> > &dist) const {
    const Real lo = dist->lowerBound(), hi = dist->upperBound();
    const Real half(0.5), one(1), N(nQuad);
    wts.clear(); wts.resize(nQuad);
    pts.clear(); pts.resize(nQuad);
    if ( hi >= one ) {
      wts[0] = half/(N-half);
      pts[0] = lo;
      for (int i = 1; i < nQuad; ++i) {
        wts[i] = one/(N-half);
        pts[i] = dist->invertCDF(static_cast<Real>(i)/N);
      }
    }
    else {
      wts[0] = half/(N-one);
      pts[0] = lo;
      for (int i = 1; i < nQuad-1; ++i) {
        wts[i] = one/(N-one);
        pts[i] = dist->invertCDF(static_cast<Real>(i)/N);
      }
      wts[nQuad-1] = half/(N-one);
      pts[nQuad-1] = hi;
    }
  }

  void printQuad(const std::vector<Real> &pts,
                 const std::vector<Real> &wts,
                 const bool print = false) const {
    if ( print ) {
      const int nQuad = wts.size();
      std::cout << std::endl;
      std::cout << std::scientific << std::setprecision(15);
      std::cout << std::setw(25) << std::left << "Points"
                << std::setw(25) << std::left << "Weights"
                << std::endl;
      for (int i = 0; i < nQuad; ++i) {
        std::cout << std::setw(25) << std::left << pts[i]
                  << std::setw(25) << std::left << wts[i]
                  << std::endl;
      }
      std::cout << std::endl;
    }
  }


public:
  SpectralRisk(void) : RandVarFunctional<Real>() {}

  SpectralRisk( const ROL::Ptr<Distribution<Real> > &dist,
                const int nQuad,
                const ROL::Ptr<PlusFunction<Real> > &pf)
    : RandVarFunctional<Real>() {
    // Build generalized trapezoidal rule
    std::vector<Real> wts(nQuad), pts(nQuad);
    buildQuadFromDist(pts,wts,nQuad,dist);
    // Build mixed quantile quadrangle risk measure
    buildMixedQuantile(pts,wts,pf);
    // Check inputs
    checkInputs(dist);
  }

  SpectralRisk(ROL::ParameterList &parlist)
    : RandVarFunctional<Real>() {
    // Parse parameter list
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Spectral Risk");
    int nQuad  = list.get("Number of Quadrature Points",5);
    bool print = list.get("Print Quadrature to Screen",false);
    // Build distribution
    ROL::Ptr<Distribution<Real> > dist = DistributionFactory<Real>(list);
    // Build plus function approximation
    ROL::Ptr<PlusFunction<Real> > pf = ROL::makePtr<PlusFunction<Real>>(list);
    // Build generalized trapezoidal rule
    std::vector<Real> wts(nQuad), pts(nQuad);
    buildQuadFromDist(pts,wts,nQuad,dist);
    printQuad(pts,wts,print);
    // Build mixed quantile quadrangle risk measure
    buildMixedQuantile(pts,wts,pf);
    // Check inputs
    checkInputs(dist);
  }

  SpectralRisk( const std::vector<Real> &pts, const std::vector<Real> &wts,
                const ROL::Ptr<PlusFunction<Real> > &pf)
    : RandVarFunctional<Real>() {
    buildMixedQuantile(pts,wts,pf);
    // Check inputs
    checkInputs();
  }

  void setStorage(const Ptr<ScalarController<Real>> &value_storage,
                  const Ptr<VectorController<Real>> &gradient_storage) {
    RandVarFunctional<Real>::setStorage(value_storage,gradient_storage);
    mqq_->setStorage(value_storage,gradient_storage);
  }

  void setHessVecStorage(const Ptr<ScalarController<Real>> &gradvec_storage,
                         const Ptr<VectorController<Real>> &hessvec_storage) {
    RandVarFunctional<Real>::setHessVecStorage(gradvec_storage,hessvec_storage);
    mqq_->setHessVecStorage(gradvec_storage,hessvec_storage);
  }

  void setSample(const std::vector<Real> &point, const Real weight) {
    RandVarFunctional<Real>::setSample(point,weight);
    mqq_->setSample(point,weight);
  }

  void resetStorage(bool flag = true) {
    RandVarFunctional<Real>::resetStorage(flag);
    mqq_->resetStorage(flag);
  }

  void resetStorage(UpdateType type) {
    RandVarFunctional<Real>::resetStorage(type);
    mqq_->resetStorage(type);
  }

  void initialize(const Vector<Real> &x) {
    RandVarFunctional<Real>::initialize(x);
    mqq_->initialize(x);
  }

  Real computeStatistic(const Ptr<const std::vector<Real>> &xstat) const override {
    return mqq_->computeStatistic(xstat);
  }

  void updateValue(Objective<Real>         &obj,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   Real                    &tol) {
    mqq_->updateValue(obj,x,xstat,tol);
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    mqq_->updateGradient(obj,x,xstat,tol);
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    mqq_->updateHessVec(obj,v,vstat,x,xstat,tol);
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    return mqq_->getValue(x,xstat,sampler);
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    mqq_->getGradient(g,gstat,x,xstat,sampler);
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    mqq_->getHessVec(hv,hvstat,v,vstat,x,xstat,sampler);
  }
};

}

#endif
