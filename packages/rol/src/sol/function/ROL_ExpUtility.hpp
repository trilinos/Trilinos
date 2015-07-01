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

#ifndef ROL_EXPUTILITY_HPP
#define ROL_EXPUTILITY_HPP

#include "ROL_RiskMeasure.hpp"

namespace ROL {

template<class Real>
class ExpUtility : public RiskMeasure<Real> {
private:
  Teuchos::RCP<Vector<Real> > scaledGradient_;
  Teuchos::RCP<Vector<Real> > dualVector_;
  bool firstReset_;

public:
  ExpUtility(void) : RiskMeasure<Real>(), firstReset_(true) {}

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    RiskMeasure<Real>::reset(x0,x);
    if ( firstReset_ ) {
      scaledGradient_ = (x.dual()).clone();
      dualVector_ = (x.dual()).clone();
      firstReset_ = false;
    }
    scaledGradient_->zero();
    dualVector_->zero();
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x,
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    RiskMeasure<Real>::reset(x0,x,v0,v);
    if ( firstReset_ ) {
      scaledGradient_ = (x.dual()).clone();
      dualVector_ = (x.dual()).clone();
      firstReset_ = false;
    }
    scaledGradient_->zero();
    dualVector_->zero();
  }

  void update(const Real val, const Real weight) {
    RiskMeasure<Real>::val_ += weight * std::exp(val);
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    Real ev = std::exp(val);
    RiskMeasure<Real>::val_ += weight * ev;
    RiskMeasure<Real>::g_->axpy(weight*ev,g);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
                      const Real weight) {
    Real ev = std::exp(val);
    RiskMeasure<Real>::val_ += weight * ev;
    RiskMeasure<Real>::gv_  -= weight * ev * gv;
    RiskMeasure<Real>::g_->axpy(weight*ev,g);
    RiskMeasure<Real>::hv_->axpy(weight*ev,hv);
    scaledGradient_->axpy(weight*ev*gv,g);
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    Real val = RiskMeasure<Real>::val_;
    Real ev  = 0.0;
    sampler.sumAll(&val,&ev,1);
    return std::log(ev);
  }

  virtual void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    Real val = RiskMeasure<Real>::val_;
    Real ev  = 0.0;
    sampler.sumAll(&val,&ev,1);
    sampler.sumAll(*(RiskMeasure<Real>::g_),g);
    g.scale(1.0/ev);
  }

  virtual void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    Real val = RiskMeasure<Real>::val_;
    Real ev  = 0.0;
    sampler.sumAll(&val,&ev,1);

    Real gv  = RiskMeasure<Real>::gv_;
    Real egv = 0.0;
    sampler.sumAll(&gv,&egv,1);

    sampler.sumAll(*(RiskMeasure<Real>::hv_),hv);

    sampler.sumAll(*scaledGradient_,*dualVector_);
    hv.plus(*dualVector_);
    hv.scale(1.0/ev);

    dualVector_->zero();
    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector_);
    dualVector_->scale(egv/(ev*ev));
    hv.plus(*dualVector_);
  }
};

}

#endif
