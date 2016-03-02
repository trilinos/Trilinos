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

#ifndef ROL_BPOEOBJECTIVE_HPP
#define ROL_BPOEOBJECTIVE_HPP

#include "Teuchos_RCP.hpp"
#include "ROL_RiskVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_ParametrizedObjective.hpp"
#include "ROL_SampleGenerator.hpp"

namespace ROL {

template<class Real>
class BPOEObjective : public Objective<Real> {
private:
  Teuchos::RCP<ParametrizedObjective<Real> > ParametrizedObjective_;

  Real order_;
  Real threshold_;

  Teuchos::RCP<SampleGenerator<Real> > ValueSampler_;
  Teuchos::RCP<SampleGenerator<Real> > GradientSampler_;
  Teuchos::RCP<SampleGenerator<Real> > HessianSampler_;

  Teuchos::RCP<Vector<Real> > pointGrad_;
  Teuchos::RCP<Vector<Real> > pointHess_;

  Teuchos::RCP<Vector<Real> > gradient0_;
  Teuchos::RCP<Vector<Real> > sumGrad0_;
  Teuchos::RCP<Vector<Real> > gradient1_;
  Teuchos::RCP<Vector<Real> > sumGrad1_;
  Teuchos::RCP<Vector<Real> > gradient2_;
  Teuchos::RCP<Vector<Real> > sumGrad2_;
  Teuchos::RCP<Vector<Real> > hessvec_;
  Teuchos::RCP<Vector<Real> > sumHess_;
 
  bool storage_;
  bool initialized_;

  std::map<std::vector<Real>,Real> value_storage_;
  std::map<std::vector<Real>,Teuchos::RCP<Vector<Real> > > gradient_storage_;

  void initialize(const Vector<Real> &x) {
    pointGrad_ = x.dual().clone();
    pointHess_ = x.dual().clone();
    gradient0_ = x.dual().clone();
    sumGrad0_  = x.dual().clone();
    gradient1_ = x.dual().clone();
    sumGrad1_  = x.dual().clone();
    gradient2_ = x.dual().clone();
    sumGrad2_  = x.dual().clone();
    hessvec_   = x.dual().clone();
    sumHess_   = x.dual().clone();
    initialized_ = true;
  }

  void unwrap_const_CVaR_vector(Teuchos::RCP<Vector<Real> > &xvec, Real &xvar,
                          const Vector<Real> &x) {
    xvec = Teuchos::rcp_const_cast<Vector<Real> >(
      Teuchos::dyn_cast<const RiskVector<Real> >(
        Teuchos::dyn_cast<const Vector<Real> >(x)).getVector());
    xvar = Teuchos::dyn_cast<const RiskVector<Real> >(
        Teuchos::dyn_cast<const Vector<Real> >(x)).getStatistic();
    if ( !initialized_ ) {
      initialize(*xvec);
    }
  }

  void getValue(Real &val, const Vector<Real> &x,
          const std::vector<Real> &param, Real &tol) {
    if ( storage_ && value_storage_.count(param) ) {
      val = value_storage_[param];
    }
    else {
      ParametrizedObjective_->setParameter(param);
      val = ParametrizedObjective_->value(x,tol);
      if ( storage_ ) {
        value_storage_.insert(std::pair<std::vector<Real>,Real>(param,val));
      }
    }
  }

  void getGradient(Vector<Real> &g, const Vector<Real> &x,
             const std::vector<Real> &param, Real &tol) {
    if ( storage_ && gradient_storage_.count(param) ) {
      g.set(*(gradient_storage_[param]));
    }
    else {
      ParametrizedObjective_->setParameter(param);
      ParametrizedObjective_->gradient(g,x,tol);
      if ( storage_ ) {
        Teuchos::RCP<Vector<Real> > tmp = g.clone();
        gradient_storage_.insert(std::pair<std::vector<Real>,Teuchos::RCP<Vector<Real> > >(param,tmp));
        gradient_storage_[param]->set(g);
      }
    }
  }

  void getHessVec(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x,
            const std::vector<Real> &param, Real &tol) {
    ParametrizedObjective_->setParameter(param);
    ParametrizedObjective_->hessVec(hv,v,x,tol);
  }

public:
  virtual ~BPOEObjective() {}

  BPOEObjective( const Teuchos::RCP<ParametrizedObjective<Real> > &pObj,
                 const Real order, const Real threshold,
                 const Teuchos::RCP<SampleGenerator<Real> >       &vsampler, 
                 const Teuchos::RCP<SampleGenerator<Real> >       &gsampler,
                 const Teuchos::RCP<SampleGenerator<Real> >       &hsampler,
                 const bool storage = true )
    : ParametrizedObjective_(pObj), order_(order), threshold_(threshold),
      ValueSampler_(vsampler), GradientSampler_(gsampler), HessianSampler_(hsampler),
      storage_(storage), initialized_(false) {
    value_storage_.clear();
    gradient_storage_.clear();
  }

  BPOEObjective( const Teuchos::RCP<ParametrizedObjective<Real> > &pObj,
                 const Real order, const Real threshold,
                 const Teuchos::RCP<SampleGenerator<Real> >       &vsampler, 
                 const Teuchos::RCP<SampleGenerator<Real> >       &gsampler,
                 const bool storage = true )
    : ParametrizedObjective_(pObj), order_(order), threshold_(threshold),
      ValueSampler_(vsampler), GradientSampler_(gsampler), HessianSampler_(gsampler),
      storage_(storage), initialized_(false) {
    value_storage_.clear();
    gradient_storage_.clear();
  }

  BPOEObjective( const Teuchos::RCP<ParametrizedObjective<Real> > &pObj,
                 const Real order, const Real threshold,
                 const Teuchos::RCP<SampleGenerator<Real> >       &sampler,
                 const bool storage = true )
    : ParametrizedObjective_(pObj), order_(order), threshold_(threshold),
      ValueSampler_(sampler), GradientSampler_(sampler), HessianSampler_(sampler),
      storage_(storage), initialized_(false) {
    value_storage_.clear();
    gradient_storage_.clear();
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    Teuchos::RCP<Vector<Real> > xvec; Real xvar = 0.0;
    unwrap_const_CVaR_vector(xvec,xvar,x);
    ParametrizedObjective_->update(*xvec,flag,iter);
    ValueSampler_->update(*xvec);
    if ( storage_ ) {
      value_storage_.clear();
    }
    if ( flag ) {
      GradientSampler_->update(*xvec);
      HessianSampler_->update(*xvec);
      if ( storage_ ) {
        gradient_storage_.clear();
      }
    }
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    Teuchos::RCP<Vector<Real> > xvec; Real xvar = 0.0;
    unwrap_const_CVaR_vector(xvec,xvar,x);
    // Initialize storage
    std::vector<Real> point;
    Real weight = 0.0, myval = 0.0, pval = 0.0, val = 0.0, bp = 0.0;
    int start = ValueSampler_->start(), end = ValueSampler_->numMySamples();
    for ( int i = start; i < end; i++ ) {
      weight = ValueSampler_->getMyWeight(i);
      point  = ValueSampler_->getMyPoint(i);
      // Compute f(xvec,xi)
      getValue(pval,*xvec,point,tol);
      bp = xvar*(pval-threshold_)+1.0;
      if ( bp > 0.0 ) {
        // Build partial sum depending on value
        myval += weight*((order_==1.0) ? bp
                          : std::pow(bp,order_));
      }
    }
    // Update expected value
    ValueSampler_->sumAll(&myval,&val,1);
    // Return BPOE value
    return ((order_==1.0) ? val : std::pow(val,1.0/order_));
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    Teuchos::RCP<Vector<Real> > xvec; Real xvar = 0.0;
    unwrap_const_CVaR_vector(xvec,xvar,x);
    RiskVector<Real> &gc = Teuchos::dyn_cast<RiskVector<Real> >(
      Teuchos::dyn_cast<Vector<Real> >(g));
    // Initialize storage
    g.zero(); sumGrad0_->zero(); pointGrad_->zero();
    std::vector<Real> point, val(2,0.0), myval(2,0.0);
    Real weight = 0.0, pval = 0.0, pvalp0 = 0.0, pvalp1 = 0.0, bp = 0.0;
    int start = GradientSampler_->start(), end = GradientSampler_->numMySamples();
    for ( int i = start; i < end; i++ ) {
      weight = GradientSampler_->getMyWeight(i);
      point  = GradientSampler_->getMyPoint(i);
      // Compute the value of f(xvec,xi)
      getValue(pval,*xvec,point,tol);
      bp = xvar*(pval-threshold_)+1.0;
      if ( bp > 0.0 ) {
        // Compute max(0,f(xvec,xi)-xvar)^order
        pvalp0 = ((order_==1.0) ? bp
                   : std::pow(bp,order_));
        pvalp1 = ((order_==1.0) ? 1.0
                   : ((order_==2.0) ? bp
                     : std::pow(bp,order_-1.0)));
        // Build partial sums depending on value
        myval[0] += weight*pvalp0;
        myval[1] += weight*pvalp1*(pval-threshold_);
        // Compute gradient of f(xvec,xi)
        getGradient(*pointGrad_,*xvec,point,tol);
        // Build partial sum depending on gradient
        sumGrad0_->axpy(weight*pvalp1,*pointGrad_);
      }
    }
    // Combine partial sums
    GradientSampler_->sumAll(&myval[0],&val[0],2);
    // Compute VaR gradient and BPOE gradient
    Real gvar = 0.0; gradient0_->zero();
    if ( std::abs(val[0]) >= ROL_EPSILON<Real>()) {
      GradientSampler_->sumAll(*sumGrad0_,*gradient0_);
      Real norm = std::pow(val[0],(order_-1.0)/order_); 
      gradient0_->scale(xvar/norm);
      gvar = val[1]/norm;
    }
    // Set gradient components of CVaR vector
    gc.setStatistic(gvar);
    gc.setVector(*(Teuchos::rcp_dynamic_cast<Vector<Real> >(gradient0_)));
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v,
                        const Vector<Real> &x, Real &tol ) {
    Teuchos::RCP<Vector<Real> > xvec; Real xvar = 0.0;
    unwrap_const_CVaR_vector(xvec,xvar,x);
    Teuchos::RCP<Vector<Real> > vvec; Real vvar = 0.0;
    unwrap_const_CVaR_vector(vvec,vvar,v);
    RiskVector<Real> &hvc = Teuchos::dyn_cast<RiskVector<Real> >(
      Teuchos::dyn_cast<Vector<Real> >(hv));
    // Initialize storage
    hv.zero(); sumHess_->zero(); hessvec_->zero();
    sumGrad0_->zero(); sumGrad1_->zero(); sumGrad2_->zero();
    gradient0_->zero(); gradient1_->zero(); gradient2_->zero();
    std::vector<Real> point, val(5,0.0), myval(5,0.0);
    Real pval = 0.0, pvalp0 = 0.0, pvalp1 = 0.0, pvalp2 = 0.0;
    Real weight = 0.0, gv = 0.0, bp = 0.0;
    int start = HessianSampler_->start(), end = HessianSampler_->numMySamples();
    for ( int i = start; i < end; i++ ) {
      // Get sample and associated probability
      weight = HessianSampler_->getMyWeight(i);
      point  = HessianSampler_->getMyPoint(i);
      // Compute the value of f(xvec,xi)
      getValue(pval,*xvec,point,tol);
      bp = xvar*(pval-threshold_)+1.0;
      if ( bp > 0.0 ) {
        // Compute max(0,f(xvec,xi)-xvar)^order
        pvalp0 = ((order_==1.0) ? bp
                   : std::pow(bp,order_));
        pvalp1 = ((order_==1.0) ? 1.0
                   : ((order_==2.0) ? bp
                     : std::pow(bp,order_-1.0)));
        pvalp2 = ((order_==1.0) ? 0.0
                   : ((order_==2.0) ? 1.0
                     : ((order_==3.0) ? bp
                       : std::pow(bp,order_-2.0))));
        // Build partial sums depending on value
        myval[0] += weight*pvalp0;
        myval[1] += weight*pvalp1*(pval-threshold_);
        myval[2] += weight*pvalp2*(pval-threshold_)*(pval-threshold_);
        // Compute the gradient and directional derivative of f(xvec,xi)
        getGradient(*pointGrad_,*xvec,point,tol);
        gv = pointGrad_->dot(vvec->dual());
        // Build partial sums depending on gradient
        myval[3] += weight*pvalp1*gv;
        myval[4] += weight*pvalp2*(pval-threshold_)*gv;
        sumGrad0_->axpy(weight*pvalp1,*pointGrad_);
        sumGrad1_->axpy(weight*pvalp2*(pval-threshold_),*pointGrad_);
        sumGrad2_->axpy(weight*pvalp2*gv,*pointGrad_);
        // Compute the hessian of f(xvec,xi) in the direction vvec
        getHessVec(*pointHess_,*vvec,*xvec,point,tol);
        // Build partial sum depending on the hessian
        sumHess_->axpy(weight*pvalp1,*pointHess_);
      }
    }
    // Compile partial sums
    HessianSampler_->sumAll(&myval[0],&val[0],5);
    Real hvar = 0.0; hessvec_->zero();
    if ( std::abs(val[0]) >= ROL_EPSILON<Real>() ) {
      HessianSampler_->sumAll(*sumGrad0_,*gradient0_);
      HessianSampler_->sumAll(*sumGrad1_,*gradient1_);
      HessianSampler_->sumAll(*sumGrad2_,*gradient2_);
      HessianSampler_->sumAll(*sumHess_,*hessvec_);
      // Compute VaR Hessian-times-a-vector and BPOE Hessian-times-a-vector
      Real norm0 = ((order_==1.0) ? 1.0
                     : ((order_==2.0) ? std::sqrt(val[0])
                       : std::pow(val[0],(order_-1.0)/order_))); 
      Real norm1 = ((order_==1.0) ? val[0]
                     : std::pow(val[0],(2.0*order_-1.0)/order_)); 
      hvar = (order_-1.0)*((val[2]/norm0 - val[1]*val[1]/norm1)*vvar
                          +xvar*(val[4]/norm0 - val[3]*val[1]/norm1))
                          +(val[3]/norm0);
      hessvec_->scale(xvar/norm0); //(order_-1.0)/norm0);
      Real coeff = -(order_-1.0)*xvar*(xvar*val[3]+vvar*val[1])/norm1+vvar/norm0;
      hessvec_->axpy(coeff,*gradient0_);
      hessvec_->axpy((order_-1.0)*vvar*xvar/norm0,*gradient1_);
      hessvec_->axpy((order_-1.0)*xvar*xvar/norm0,*gradient2_);
    }
    // Set gradient components of CVaR vector
    hvc.setStatistic(hvar);
    hvc.setVector(*(Teuchos::rcp_dynamic_cast<Vector<Real> >(hessvec_)));
  }

  virtual void precond( Vector<Real> &Pv, const Vector<Real> &v,
                        const Vector<Real> &x, Real &tol ) {
    Pv.set(v.dual());
  }
};

}

#endif
