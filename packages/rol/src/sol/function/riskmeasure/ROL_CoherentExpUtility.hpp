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

#ifndef ROL_COHERENTEXPUTILITY_HPP
#define ROL_COHERENTEXPUTILITY_HPP

#include "ROL_RiskMeasure.hpp"

/** @ingroup risk_group
    \class ROL::CoherentExpUtility
    \brief Provides the interface for the coherent entropic risk measure.

    The coherent entropic risk measure is
    \f[
       \mathcal{R}(X) = \inf_{\lambda > 0} \left\{
         \lambda \log\mathbb{E}\left[\exp\left(\frac{X}{\lambda}\right)\right]
         \right\}.
    \f]
    \f$\mathcal{R}\f$ is a law-invariant coherent risk measure.
    ROL implements this by augmenting the optimization vector \f$x_0\f$ with
    the parameter \f$\lambda\f$, then minimizes jointly for \f$(x_0,\lambda)\f$.
*/

namespace ROL {

template<class Real>
class CoherentExpUtility : public RiskMeasure<Real> {
private:
  bool firstReset_;

  Teuchos::RCP<Vector<Real> > scaledGradient1_;
  Teuchos::RCP<Vector<Real> > scaledGradient2_;
  Real dval1_;
  Real dval2_;
  Real dval3_;

  Teuchos::RCP<Vector<Real> > dualVector1_;
  Teuchos::RCP<Vector<Real> > dualVector2_;

  Real xstat_;
  Real vstat_;

public:
  CoherentExpUtility(void) : RiskMeasure<Real>(), firstReset_(true),
    dval1_(0), dval2_(0), dval3_(0), xstat_(0), vstat_(0) {}

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    Real zero(0);
    RiskMeasure<Real>::reset(x0,x);
    xstat_ = Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic();
    if ( firstReset_ ) {
      scaledGradient1_ = (x0->dual()).clone();
      scaledGradient2_ = (x0->dual()).clone();
      dualVector1_ = (x0->dual()).clone();
      dualVector2_ = (x0->dual()).clone();
      firstReset_ = false;
    }
    scaledGradient1_->zero(); scaledGradient2_->zero();
    dualVector1_->zero(); dualVector2_->zero();
    dval1_ = zero; dval2_ = zero; dval3_ = zero;
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x,
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    reset(x0,x);
    v0 = Teuchos::rcp_const_cast<Vector<Real> >(
           Teuchos::dyn_cast<const RiskVector<Real> >(v).getVector());
    vstat_ = Teuchos::dyn_cast<const RiskVector<Real> >(v).getStatistic();
  }

  void update(const Real val, const Real weight) {
    RiskMeasure<Real>::val_ += weight * std::exp(val/xstat_);
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    Real val = RiskMeasure<Real>::val_, ev(0);
    sampler.sumAll(&val,&ev,1);
    return xstat_*std::log(ev);
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    Real ev = std::exp(val/xstat_);
    RiskMeasure<Real>::val_ += weight * ev;
    RiskMeasure<Real>::gv_  += weight * ev * val;
    RiskMeasure<Real>::g_->axpy(weight*ev,g);
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    Real one(1);
    // Perform sum over batches
    std::vector<Real> myval(2,0), val(2,0);
    myval[0] = RiskMeasure<Real>::val_;
    myval[1] = RiskMeasure<Real>::gv_;
    sampler.sumAll(&myval[0],&val[0],2);
    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector1_);
    // Compute partial derivatives
    Real gstat = std::log(myval[0]) - myval[1]/(myval[0]*xstat_);
    dualVector1_->scale(one/myval[0]);
    // Set partial derivatives in g vector
    (Teuchos::dyn_cast<RiskVector<Real> >(g)).setVector(*dualVector1_);
    (Teuchos::dyn_cast<RiskVector<Real> >(g)).setStatistic(gstat);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
                      const Real weight) {
    Real ev = std::exp(val/xstat_);
    RiskMeasure<Real>::val_ += weight * ev;
    RiskMeasure<Real>::gv_  += weight * ev * gv;
    dval1_                  += weight * ev * val;
    dval2_                  += weight * ev * val * val;
    dval3_                  += weight * ev * val * gv;
    RiskMeasure<Real>::g_->axpy(weight*ev,g);
    RiskMeasure<Real>::hv_->axpy(weight*ev,hv);
    scaledGradient1_->axpy(weight*ev*gv,g);
    scaledGradient2_->axpy(weight*ev*val,g);
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    Real one(1);
    std::vector<Real> myval(5,0), val(5,0);
    myval[0] = RiskMeasure<Real>::val_;
    myval[1] = RiskMeasure<Real>::gv_;
    myval[2] = dval1_;
    myval[3] = dval2_;
    myval[4] = dval3_;
    sampler.sumAll(&myval[0],&val[0],5);

    Real xs2 = xstat_*xstat_;
    Real xs3 = xs2*xstat_;
    Real v02 = val[0]*val[0];
    Real h11 = (val[3]*val[0] - val[2]*val[2])/(v02*xs3) * vstat_;
    Real h12 = (val[1]*val[2] - val[4]*val[0])/(v02*xs2); 

    sampler.sumAll(*(RiskMeasure<Real>::hv_),*dualVector1_);
    sampler.sumAll(*scaledGradient1_,*dualVector2_);
    dualVector1_->axpy(one/xstat_,*dualVector2_);
    dualVector1_->scale(one/val[0]);
    dualVector2_->zero();
    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector2_);
    dualVector1_->axpy(vstat_*val[2]/(xs2*v02)-val[1]/(v02*xstat_),*dualVector2_);
    dualVector2_->zero();
    sampler.sumAll(*scaledGradient2_,*dualVector2_);
    dualVector1_->axpy(-vstat_/val[0],*dualVector2_);

    (Teuchos::dyn_cast<RiskVector<Real> >(hv)).setVector(*dualVector1_);
    (Teuchos::dyn_cast<RiskVector<Real> >(hv)).setStatistic(h11+h12);
  }
};

}

#endif
