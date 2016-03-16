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

namespace ROL {

template<class Real>
class KLDivergence : public RiskMeasure<Real> {
private:
  Real eps_;

  Real gval_;
  Real gvval_;
  Real hval_;
  Teuchos::RCP<Vector<Real> > scaledGradient_;
  Teuchos::RCP<Vector<Real> > scaledHessVec_;
  Teuchos::RCP<Vector<Real> > dualVector1_;
  Teuchos::RCP<Vector<Real> > dualVector2_;

  Real xstat_;
  Real vstat_;

  bool firstReset_;

public:
  KLDivergence(const Real eps = 1.e-2)
    : RiskMeasure<Real>(), firstReset_(true) {
    Real zero(0), oem2(1.e-2);
    eps_ = eps > zero ? eps : oem2;
  }

  KLDivergence(Teuchos::ParameterList &parlist)
    : RiskMeasure<Real>(), firstReset_(true) {
    Real zero(0), oem2(1.e-2);
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("KL Divergence");
    Real eps = list.get("Threshold",oem2);
    eps_ = eps > zero ? eps : oem2;
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    Real zero(0);
    RiskMeasure<Real>::reset(x0,x);
    xstat_ = Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic();
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

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x,
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    reset(x0,x);
    v0 = Teuchos::rcp_const_cast<Vector<Real> >(Teuchos::dyn_cast<const RiskVector<Real> >(v).getVector());
    vstat_ = Teuchos::dyn_cast<const RiskVector<Real> >(v).getStatistic();
  }

  void update(const Real val, const Real weight) {
    RiskMeasure<Real>::val_ += weight * std::exp(val/xstat_);
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    Real ev = std::exp(val/xstat_);
    RiskMeasure<Real>::val_ += weight * ev;
    RiskMeasure<Real>::g_->axpy(weight*ev,g);
    gval_ += weight*ev*val;
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
                      const Real weight) {
    Real ev = std::exp(val/xstat_);
    RiskMeasure<Real>::val_ += weight * ev;
    RiskMeasure<Real>::gv_  += weight * ev * gv;
    RiskMeasure<Real>::g_->axpy(weight*ev,g);
    RiskMeasure<Real>::hv_->axpy(weight*ev,hv);
    scaledGradient_->axpy(weight*ev*gv,g);
    scaledHessVec_->axpy(weight*ev*val,g);
    gval_  += weight*ev*val;
    gvval_ += weight*ev*val*gv;
    hval_  += weight*ev*val*val;
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    Real val = RiskMeasure<Real>::val_, ev(0);
    sampler.sumAll(&val,&ev,1);
    return xstat_*(eps_ + std::log(ev));
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    Real one(1);
    std::vector<Real> local(2), global(2);
    local[0] = RiskMeasure<Real>::val_;
    local[1] = gval_;
    sampler.sumAll(&local[0],&global[0],2);
    Real ev    = global[0];
    Real egval = global[1];

    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector1_);
    dualVector1_->scale(one/ev);

    Real gstat = eps_ + std::log(ev) - egval/(ev*xstat_);

    (Teuchos::dyn_cast<RiskVector<Real> >(g)).setVector(*dualVector1_);
    (Teuchos::dyn_cast<RiskVector<Real> >(g)).setStatistic(gstat);
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    Real one(1);
    std::vector<Real> local(5), global(5);
    local[0] = RiskMeasure<Real>::val_;
    local[1] = RiskMeasure<Real>::gv_;
    local[2] = gval_;
    local[3] = gvval_;
    local[4] = hval_;
    sampler.sumAll(&local[0],&global[0],5);
    Real ev     = global[0];
    Real egv    = global[1];
    Real egval  = global[2];
    Real egvval = global[3];
    Real ehval  = global[4];

    sampler.sumAll(*(RiskMeasure<Real>::hv_),*dualVector1_);

    sampler.sumAll(*scaledGradient_,*dualVector2_);
    dualVector1_->axpy(one/xstat_,*dualVector2_);
    dualVector1_->scale(one/ev);

    dualVector2_->zero();
    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector2_);
    dualVector1_->axpy((vstat_*egval/xstat_ - egv)/(xstat_*ev*ev),*dualVector2_);

    dualVector2_->zero();
    sampler.sumAll(*scaledHessVec_,*dualVector2_);
    dualVector1_->axpy(-vstat_/(xstat_*xstat_*ev),*dualVector2_);

    (Teuchos::dyn_cast<RiskVector<Real> >(hv)).setVector(*dualVector1_);

    Real hstat = vstat_/(xstat_*xstat_*xstat_*ev) * (ehval - egval*egval/ev)
                 + (egv*egval/ev - egvval)/(ev*xstat_*xstat_);
    (Teuchos::dyn_cast<RiskVector<Real> >(hv)).setStatistic(hstat);
  }
};

}

#endif
