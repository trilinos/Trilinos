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

#ifndef ROL_KLDIVERGENCE_HPP
#define ROL_KLDIVERGENCE_HPP

#include "ROL_RiskMeasure.hpp"

/** @ingroup risk_group
    \class ROL::KLDivergence
    \brief Provides an interface for the Kullback-Leibler distributionally robust
    expectation.

    This class defines a risk measure \f$\mathcal{R}\f$ which arises in distributionally
    robust stochastic programming.  \f$\mathcal{R}\f$ is given by
    \f[
       \mathcal{R}(X) = \sup_{\vartheta\in\mathfrak{A}}
           \mathbb{E}[\vartheta X]
    \f]
    where \f$\mathfrak{A}\f$ is called the ambiguity (or uncertainty) set and 
    is defined by a constraint on the Kullback-Leibler divergence, i.e.,
    \f[
       \mathfrak{A} = \{\vartheta\in\mathcal{X}^*\,:\,
         \mathbb{E}[\vartheta] = 1,\; \vartheta \ge 0,\;\text{and}\;
         \mathbb{E}[\vartheta\log(\vartheta)] \le \epsilon\}.
    \f]
    \f$\mathcal{R}\f$ is a law-invariant, coherent risk measure.  Moreover, by a
    duality argument, \f$\mathcal{R}\f$ can be reformulated as
    \f[
       \mathcal{R}(X) = \inf_{\lambda > 0}\left\{
             \lambda \epsilon + \lambda\mathbb{E}\left[\exp\left(
                \frac{X}{\lambda}\right)\right]\right\}.
    \f]
    ROL implements this by augmenting the optimization vector \f$x_0\f$ with
    the parameter \f$\lambda\f$, then minimizes jointly for \f$(x_0,\lambda)\f$.
*/

namespace ROL {

template<class Real>
class KLDivergence : public RiskMeasure<Real> {
private:
  Real eps_;

  Real gval_;
  Real gvval_;
  Real hval_;
  ROL::SharedPointer<Vector<Real> > scaledGradient_;
  ROL::SharedPointer<Vector<Real> > scaledHessVec_;
  ROL::SharedPointer<Vector<Real> > dualVector1_;
  ROL::SharedPointer<Vector<Real> > dualVector2_;

  Real xstat_;
  Real vstat_;

  bool firstReset_;

  void checkInputs(void) const {
    Real zero(0);
    TEUCHOS_TEST_FOR_EXCEPTION((eps_ <= zero), std::invalid_argument,
      ">>> ERROR (ROL::KLDivergence): Threshold must be positive!");
  }

public:
  /** \brief Constructor.

      @param[in]     eps    is the tolerance for the KL divergence constraint
  */
  KLDivergence(const Real eps = 1.e-2)
    : RiskMeasure<Real>(), eps_(eps), firstReset_(true) {
    checkInputs();
  }

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"KL Divergence" and
      within the "KL Divergence" sublist should have the following parameters
      \li "Threshold" (greater than 0)
  */
  KLDivergence(ROL::ParameterList &parlist)
    : RiskMeasure<Real>(), firstReset_(true) {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("KL Divergence");
    eps_ = list.get<Real>("Threshold");
    checkInputs();
  }

  void reset(ROL::SharedPointer<Vector<Real> > &x0, const Vector<Real> &x) {
    Real zero(0);
    RiskMeasure<Real>::reset(x0,x);
    int index = RiskMeasure<Real>::getIndex();
    int comp  = RiskMeasure<Real>::getComponent();
    xstat_ = (*dynamic_cast<const RiskVector<Real>&>(x).getStatistic(comp,index))[0];
    if ( firstReset_ ) {
      scaledGradient_ = (x0->dual()).clone();
      scaledHessVec_  = (x0->dual()).clone();
      dualVector1_ = (x0->dual()).clone();
      dualVector2_ = (x0->dual()).clone();
      firstReset_ = false;
    }
    gval_ = zero; gvval_ = zero; hval_ = zero;
    scaledGradient_->zero(); scaledHessVec_->zero();
    dualVector1_->zero(); dualVector2_->zero();
  }

  void reset(ROL::SharedPointer<Vector<Real> > &x0, const Vector<Real> &x,
             ROL::SharedPointer<Vector<Real> > &v0, const Vector<Real> &v) {
    reset(x0,x);
    v0 = ROL::constPointerCast<Vector<Real> >(dynamic_cast<const RiskVector<Real>&>(v).getVector());
    int index = RiskMeasure<Real>::getIndex();
    int comp  = RiskMeasure<Real>::getComponent();
    vstat_ = (*dynamic_cast<const RiskVector<Real>&>(v).getStatistic(comp,index))[0];
  }

  void update(const Real val, const Real weight) {
    Real ev = exponential(val,xstat_*eps_);
    RiskMeasure<Real>::val_ += weight * ev;
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    if ( xstat_ == static_cast<Real>(0) ) {
      return ROL_INF<Real>();
    }
    Real val = RiskMeasure<Real>::val_, ev(0);
    sampler.sumAll(&val,&ev,1);
    return (static_cast<Real>(1) + std::log(ev)/eps_)/xstat_;
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    Real ev = exponential(val,xstat_*eps_);
    RiskMeasure<Real>::val_ += weight * ev;
    gval_                   += weight * ev * val;
    RiskMeasure<Real>::g_->axpy(weight*ev,g);
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    std::vector<Real> local(2), global(2);
    local[0] = RiskMeasure<Real>::val_;
    local[1] = gval_;
    sampler.sumAll(&local[0],&global[0],2);
    Real ev = global[0], egval = global[1];

    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector1_);
    dualVector1_->scale(static_cast<Real>(1)/ev);
    dynamic_cast<RiskVector<Real>&>(g).setVector(*dualVector1_);

    Real gstat(0);
    if ( xstat_ == static_cast<Real>(0) ) {
      gstat = ROL_INF<Real>();
    }
    else {
      gstat = -((static_cast<Real>(1) + std::log(ev)/eps_)/xstat_
              - egval/ev)/xstat_;
    }
    int index = RiskMeasure<Real>::getIndex();
    int comp  = RiskMeasure<Real>::getComponent();
    dynamic_cast<RiskVector<Real>&>(g).setStatistic(gstat,comp,index);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
                      const Real weight) {
    Real ev = exponential(val,xstat_*eps_);
    RiskMeasure<Real>::val_ += weight * ev;
    RiskMeasure<Real>::gv_  += weight * ev * gv;
    gval_                   += weight * ev * val;
    gvval_                  += weight * ev * val * gv;
    hval_                   += weight * ev * val * val;
    RiskMeasure<Real>::g_->axpy(weight*ev,g);
    RiskMeasure<Real>::hv_->axpy(weight*ev,hv);
    scaledGradient_->axpy(weight*ev*gv,g);
    scaledHessVec_->axpy(weight*ev*val,g);
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    std::vector<Real> local(5), global(5);
    local[0] = RiskMeasure<Real>::val_;
    local[1] = RiskMeasure<Real>::gv_;
    local[2] = gval_;
    local[3] = gvval_;
    local[4] = hval_;
    sampler.sumAll(&local[0],&global[0],5);
    Real ev     = global[0], egv   = global[1], egval = global[2];
    Real egvval = global[3], ehval = global[4];
    Real c0 = static_cast<Real>(1)/ev, c1 = c0*egval, c2 = c0*egv, c3 = eps_*c0;

    sampler.sumAll(*(RiskMeasure<Real>::hv_),*dualVector1_);
    sampler.sumAll(*scaledGradient_,*dualVector2_);
    dualVector1_->axpy(xstat_*eps_,*dualVector2_);
    dualVector1_->scale(c0);

    dualVector2_->zero();
    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector2_);
    dualVector1_->axpy(-c3*(vstat_*c1 + xstat_*c2),*dualVector2_);

    dualVector2_->zero();
    sampler.sumAll(*scaledHessVec_,*dualVector2_);
    dualVector1_->axpy(vstat_*c3,*dualVector2_);

    dynamic_cast<RiskVector<Real>&>(hv).setVector(*dualVector1_);

    Real hstat(0);
    if ( xstat_ == static_cast<Real>(0) ) {
      hstat = ROL_INF<Real>();
    }
    else {
      Real xstat2 = static_cast<Real>(2)/(xstat_*xstat_);
      Real h11 = xstat2*((static_cast<Real>(1) + std::log(ev)/eps_)/xstat_ - c1)
                 + (c3*ehval - eps_*c1*c1)/xstat_;
      hstat = vstat_ * h11 + (c3*egvval - eps_*c1*c2);
    }

    int index = RiskMeasure<Real>::getIndex();
    int comp  = RiskMeasure<Real>::getComponent();
    dynamic_cast<RiskVector<Real>&>(hv).setStatistic(hstat,comp,index);
  }

private:
  Real exponential(const Real arg1, const Real arg2) const {
    if ( arg1 < arg2 ) {
      return power(exponential(arg1),arg2);
    }
    else {
      return power(exponential(arg2),arg1);
    }
  }

  Real exponential(const Real arg) const {
    if ( arg >= std::log(ROL_INF<Real>()) ) {
      return ROL_INF<Real>();
    }
    else {
      return std::exp(arg);
    }
  }

  Real power(const Real arg, const Real pow) const {
    if ( arg >= std::pow(ROL_INF<Real>(),static_cast<Real>(1)/pow) ) {
      return ROL_INF<Real>();
    }
    else {
      return std::pow(arg,pow);
    }
  }
};

}

#endif
