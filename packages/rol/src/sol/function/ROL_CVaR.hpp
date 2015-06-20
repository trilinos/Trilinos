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

#ifndef ROL_CVAR_HPP
#define ROL_CVAR_HPP

#include "ROL_RiskMeasure.hpp"
#include "ROL_PlusFunction.hpp"
#include "ROL_CVaRVector.hpp"

namespace ROL {

template<class Real>
class CVaR : public RiskMeasure<Real> {
private:
  Real prob_;
  Real coeff_;
  Real xvar_;
  Real vvar_;
  Teuchos::RCP<PlusFunction<Real> > plusFunction_;
 
public:

  CVaR( Real prob, Real coeff, Teuchos::RCP<PlusFunction<Real> > &pf ) : plusFunction_(pf) {
    prob_ = prob;
    if ( prob < 0.0 || prob > 1.0 ) {
      prob_ = 0.5;
    }
    coeff_ = coeff;
    if ( coeff < 0.0 || coeff > 1.0 ) {
      coeff_ = 1.0;
    }
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    x0 = Teuchos::rcp_const_cast<Vector<Real> >(Teuchos::dyn_cast<const CVaRVector<Real> >(
           Teuchos::dyn_cast<const Vector<Real> >(x)).getVector());
    this->xvar_ = Teuchos::dyn_cast<const CVaRVector<Real> >(Teuchos::dyn_cast<const Vector<Real> >(x)).getVaR();
    RiskMeasure<Real>::val_ = 0.0;
    RiskMeasure<Real>::g_  = x0->clone(); RiskMeasure<Real>::g_->zero();
    RiskMeasure<Real>::hv_ = x0->clone(); RiskMeasure<Real>::hv_->zero();
  }   

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x, 
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    x0 = Teuchos::rcp_const_cast<Vector<Real> >(Teuchos::dyn_cast<const CVaRVector<Real> >(
           Teuchos::dyn_cast<const Vector<Real> >(x)).getVector());
    this->xvar_ = Teuchos::dyn_cast<const CVaRVector<Real> >(Teuchos::dyn_cast<const Vector<Real> >(x)).getVaR();
    v0 = Teuchos::rcp_const_cast<Vector<Real> >(Teuchos::dyn_cast<const CVaRVector<Real> >(
           Teuchos::dyn_cast<const Vector<Real> >(v)).getVector());
    this->vvar_ = Teuchos::dyn_cast<const CVaRVector<Real> >(Teuchos::dyn_cast<const Vector<Real> >(v)).getVaR();
    RiskMeasure<Real>::val_ = 0.0;
    RiskMeasure<Real>::g_  = x0->clone(); RiskMeasure<Real>::g_->zero();
    RiskMeasure<Real>::hv_ = x0->clone(); RiskMeasure<Real>::hv_->zero();
  }   

  void update(const Real val, const Real weight) {
    Real pf = this->plusFunction_->evaluate(val-this->xvar_,0);
    RiskMeasure<Real>::val_ += weight * ((1.0-this->coeff_) * val + this->coeff_/(1.0-this->prob_) * pf);
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    Real pf = this->plusFunction_->evaluate(val-this->xvar_,1);
    RiskMeasure<Real>::val_ += weight * pf;
    Real c  = (1.0-this->coeff_) + this->coeff_/(1.0-this->prob_)*pf;
    RiskMeasure<Real>::g_->axpy(weight*c,g);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
              const Real weight) {
    Real pf1 = this->plusFunction_->evaluate(val-this->xvar_,1);
    Real pf2 = this->plusFunction_->evaluate(val-this->xvar_,2);
    RiskMeasure<Real>::val_ += weight * pf2 * (this->vvar_ - gv);
    Real c  = pf2 * this->coeff_/(1.0-this->prob_) * (-this->vvar_ + gv);
    RiskMeasure<Real>::hv_->axpy(weight*c,g);
    c = (1.0-this->coeff_) + this->coeff_/(1.0-this->prob_) * pf1;
    RiskMeasure<Real>::hv_->axpy(weight*c,hv);
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    Real val  = RiskMeasure<Real>::val_;
    Real cvar = 0.0;
    sampler.sumAll(&val,&cvar,1);
    cvar += this->coeff_ * this->xvar_;
    return cvar;
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    CVaRVector<Real> &gs = Teuchos::dyn_cast<CVaRVector<Real> >(Teuchos::dyn_cast<Vector<Real> >(g));
    Real val = RiskMeasure<Real>::val_;
    Real var = 0.0;
    sampler.sumAll(&val,&var,1);
    
    Teuchos::RCP<Vector<Real> > gz = RiskMeasure<Real>::g_->clone();
    sampler.sumAll(*(RiskMeasure<Real>::g_),*gz);
    var *= -this->coeff_/(1.0-this->prob_);
    var += this->coeff_;
    gs.setVaR(var);
    gs.setVector(*(Teuchos::rcp_dynamic_cast<Vector<Real> >(gz))); 
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    CVaRVector<Real> &hs = Teuchos::dyn_cast<CVaRVector<Real> >(Teuchos::dyn_cast<Vector<Real> >(hv));
    Real val = RiskMeasure<Real>::val_;
    Real var = 0.0;
    sampler.sumAll(&val,&var,1);

    Teuchos::RCP<Vector<Real> > hz = RiskMeasure<Real>::hv_->clone();
    sampler.sumAll(*(RiskMeasure<Real>::hv_),*hz);
    var *= this->coeff_/(1.0-this->prob_);
    hs.setVaR(var);
    hs.setVector(*(Teuchos::rcp_dynamic_cast<Vector<Real> >(hz)));
  }
};

}

#endif
