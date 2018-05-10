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

#ifndef ROL_RISKNEUTRALOBJECTIVE_HPP
#define ROL_RISKNEUTRALOBJECTIVE_HPP

#include "Teuchos_RefCountPtr.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_SampleGenerator.hpp"

namespace ROL {

template<class Real>
class RiskNeutralObjective : public Objective<Real> {
private:
  ROL::Ptr<Objective<Real> >       ParametrizedObjective_;
  ROL::Ptr<SampleGenerator<Real> > ValueSampler_;
  ROL::Ptr<SampleGenerator<Real> > GradientSampler_;
  ROL::Ptr<SampleGenerator<Real> > HessianSampler_;

  Real value_;
  ROL::Ptr<Vector<Real> > gradient_;
  ROL::Ptr<Vector<Real> > pointDual_;
  ROL::Ptr<Vector<Real> > sumDual_;

  bool firstUpdate_;
  bool storage_;

  std::map<std::vector<Real>,Real> value_storage_;
  std::map<std::vector<Real>,ROL::Ptr<Vector<Real> > > gradient_storage_;

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
        ROL::Ptr<Vector<Real> > tmp = g.clone();
        gradient_storage_.insert(std::pair<std::vector<Real>,ROL::Ptr<Vector<Real> > >(param,tmp));
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
  virtual ~RiskNeutralObjective() {}

  RiskNeutralObjective( const ROL::Ptr<Objective<Real> >       &pObj,
                        const ROL::Ptr<SampleGenerator<Real> > &vsampler, 
                        const ROL::Ptr<SampleGenerator<Real> > &gsampler,
                        const ROL::Ptr<SampleGenerator<Real> > &hsampler,
                        const bool storage = true )
    : ParametrizedObjective_(pObj),
      ValueSampler_(vsampler), GradientSampler_(gsampler), HessianSampler_(hsampler),
      firstUpdate_(true), storage_(storage) {
    value_storage_.clear();
    gradient_storage_.clear();
  }

  RiskNeutralObjective( const ROL::Ptr<Objective<Real> >       &pObj,
                        const ROL::Ptr<SampleGenerator<Real> > &vsampler, 
                        const ROL::Ptr<SampleGenerator<Real> > &gsampler,
                        const bool storage = true )
    : ParametrizedObjective_(pObj),
      ValueSampler_(vsampler), GradientSampler_(gsampler), HessianSampler_(gsampler),
      firstUpdate_(true), storage_(storage) {
    value_storage_.clear();
    gradient_storage_.clear();
  }

  RiskNeutralObjective( const ROL::Ptr<Objective<Real> >       &pObj,
                        const ROL::Ptr<SampleGenerator<Real> > &sampler,
                        const bool storage = true )
    : ParametrizedObjective_(pObj),
      ValueSampler_(sampler), GradientSampler_(sampler), HessianSampler_(sampler),
      firstUpdate_(true), storage_(storage) {
    value_storage_.clear();
    gradient_storage_.clear();
  }

  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    if ( firstUpdate_ ) {
      gradient_  = (x.dual()).clone();
      pointDual_ = (x.dual()).clone();
      sumDual_   = (x.dual()).clone();
      firstUpdate_ = false;
    }
    ParametrizedObjective_->update(x,(flag && iter>=0),iter);
    ValueSampler_->update(x);
    value_ = static_cast<Real>(0);
    if ( storage_ ) {
      value_storage_.clear();
    }
    if ( flag && iter>=0 ) {
      GradientSampler_->update(x);
      HessianSampler_->update(x);
      gradient_->zero();
      if ( storage_ ) {
        gradient_storage_.clear();
      }
    }
  }

  virtual Real value( const Vector<Real> &x, Real &tol ) {
    Real myval(0), ptval(0), val(0), one(1), two(2), error(two*tol + one);
    std::vector<Real> ptvals;
    while ( error > tol ) {
      ValueSampler_->refine();
      for ( int i = ValueSampler_->start(); i < ValueSampler_->numMySamples(); ++i ) {
        getValue(ptval,x,ValueSampler_->getMyPoint(i),tol);
        myval += ValueSampler_->getMyWeight(i)*ptval;
        ptvals.push_back(ptval);
      }
      error = ValueSampler_->computeError(ptvals);
      ptvals.clear();
    }
    ValueSampler_->sumAll(&myval,&val,1);
    value_ += val;
    ValueSampler_->setSamples();
    tol = error;
    return value_;
  }

  virtual void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    g.zero(); pointDual_->zero(); sumDual_->zero();
    std::vector<ROL::Ptr<Vector<Real> > > ptgs;
    Real one(1), two(2), error(two*tol + one);
    while ( error > tol ) {
      GradientSampler_->refine();
      for ( int i = GradientSampler_->start(); i < GradientSampler_->numMySamples(); ++i ) {
        getGradient(*pointDual_,x,GradientSampler_->getMyPoint(i),tol);
        sumDual_->axpy(GradientSampler_->getMyWeight(i),*pointDual_);
        ptgs.push_back(pointDual_->clone());
        (ptgs.back())->set(*pointDual_);
      }
      error = GradientSampler_->computeError(ptgs,x);
//if (GradientSampler_->batchID()==0) {
//  std::cout << "IN GRADIENT: ERROR = " << error << "  TOL = " << tol << std::endl;  
//}
      ptgs.clear();
    }
    GradientSampler_->sumAll(*sumDual_,g);
    gradient_->plus(g);
    g.set(*(gradient_));
    GradientSampler_->setSamples();
    tol = error;
  }

  virtual void hessVec( Vector<Real> &hv, const Vector<Real> &v,
                        const Vector<Real> &x, Real &tol ) {
    hv.zero(); pointDual_->zero(); sumDual_->zero();
    for ( int i = 0; i < HessianSampler_->numMySamples(); ++i ) {
      getHessVec(*pointDual_,v,x,HessianSampler_->getMyPoint(i),tol);
      sumDual_->axpy(HessianSampler_->getMyWeight(i),*pointDual_);
    }
    HessianSampler_->sumAll(*sumDual_,hv);
  }

  virtual void precond( Vector<Real> &Pv, const Vector<Real> &v,
                        const Vector<Real> &x, Real &tol ) {
    Pv.set(v.dual());
  }
};

}

#endif
