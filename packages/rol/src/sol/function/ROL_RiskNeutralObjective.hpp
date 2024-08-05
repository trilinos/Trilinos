// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_RISKNEUTRALOBJECTIVE_HPP
#define ROL_RISKNEUTRALOBJECTIVE_HPP

#include "ROL_Vector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_SampleGenerator.hpp"
#include "ROL_ScalarController.hpp"
#include "ROL_VectorController.hpp"

namespace ROL {

template<class Real>
class RiskNeutralObjective : public Objective<Real> {
private:
  Ptr<Objective<Real>>       ParametrizedObjective_;
  Ptr<SampleGenerator<Real>> ValueSampler_;
  Ptr<SampleGenerator<Real>> GradientSampler_;
  Ptr<SampleGenerator<Real>> HessianSampler_;

  Real value_;
  Ptr<Vector<Real>> gradient_;
  Ptr<Vector<Real>> pointDual_;
  Ptr<Vector<Real>> sumDual_;

  bool firstUpdate_;
  bool storage_;

  //std::map<std::vector<Real>,Real> value_storage_;
  //std::map<std::vector<Real>,Ptr<Vector<Real>>> gradient_storage_;
  Ptr<ScalarController<Real>> value_storage_;
  Ptr<VectorController<Real>> gradient_storage_;

  void initialize(const Vector<Real> &x) {
    if ( firstUpdate_ ) {
      gradient_  = (x.dual()).clone();
      pointDual_ = (x.dual()).clone();
      sumDual_   = (x.dual()).clone();
      firstUpdate_ = false;
    }
  }

  void getValue(Real &val, const Vector<Real> &x,
          const std::vector<Real> &param, Real &tol) {
    bool isComputed = false;
    if ( storage_) {
      isComputed = value_storage_->get(val,param);
    }
    if (!isComputed || !storage_) {
      ParametrizedObjective_->setParameter(param);
      val = ParametrizedObjective_->value(x,tol);
      if ( storage_ ) {
        value_storage_->set(val,param);
      }
    }
  }

  void getGradient(Vector<Real> &g, const Vector<Real> &x,
             const std::vector<Real> &param, Real &tol) {
    bool isComputed = false;
    if ( storage_) {
      isComputed = gradient_storage_->get(g,param);
    }
    if (!isComputed || !storage_) {
      ParametrizedObjective_->setParameter(param);
      ParametrizedObjective_->gradient(g,x,tol);
      if ( storage_ ) {
        gradient_storage_->set(g,param);
      }
    }
  }

  void getHessVec(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x,
            const std::vector<Real> &param, Real &tol) {
    ParametrizedObjective_->setParameter(param);
    ParametrizedObjective_->hessVec(hv,v,x,tol);
  }


public:
  RiskNeutralObjective( const Ptr<Objective<Real>>       &pObj,
                        const Ptr<SampleGenerator<Real>> &vsampler, 
                        const Ptr<SampleGenerator<Real>> &gsampler,
                        const Ptr<SampleGenerator<Real>> &hsampler,
                        const bool storage = true )
    : ParametrizedObjective_(pObj),
      ValueSampler_(vsampler), GradientSampler_(gsampler), HessianSampler_(hsampler),
      firstUpdate_(true), storage_(storage) {
    value_storage_ = makePtr<ScalarController<Real>>();
    gradient_storage_ = makePtr<VectorController<Real>>();
  }

  RiskNeutralObjective( const Ptr<Objective<Real>>       &pObj,
                        const Ptr<SampleGenerator<Real>> &vsampler, 
                        const Ptr<SampleGenerator<Real>> &gsampler,
                        const bool storage = true )
    : ParametrizedObjective_(pObj),
      ValueSampler_(vsampler), GradientSampler_(gsampler), HessianSampler_(gsampler),
      firstUpdate_(true), storage_(storage) {
    value_storage_ = makePtr<ScalarController<Real>>();
    gradient_storage_ = makePtr<VectorController<Real>>();
  }

  RiskNeutralObjective( const Ptr<Objective<Real>>       &pObj,
                        const Ptr<SampleGenerator<Real>> &sampler,
                        const bool storage = true )
    : ParametrizedObjective_(pObj),
      ValueSampler_(sampler), GradientSampler_(sampler), HessianSampler_(sampler),
      firstUpdate_(true), storage_(storage) {
    value_storage_ = makePtr<ScalarController<Real>>();
    gradient_storage_ = makePtr<VectorController<Real>>();
  }

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) {
    initialize(x);
//    ParametrizedObjective_->update(x,(flag && iter>=0),iter);
    ParametrizedObjective_->update(x,type,iter);
    ValueSampler_->update(x);
    value_ = static_cast<Real>(0);
    if ( storage_ ) {
      value_storage_->objectiveUpdate(type);
      gradient_storage_->objectiveUpdate(type);
    }
    if ( type != UpdateType::Trial && type != UpdateType::Revert ) { //&& iter>=0 ) {
      GradientSampler_->update(x);
      HessianSampler_->update(x);
      gradient_->zero();
    }
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    initialize(x);
//    ParametrizedObjective_->update(x,(flag && iter>=0),iter);
    ParametrizedObjective_->update(x,flag,iter);
    ValueSampler_->update(x);
    value_ = static_cast<Real>(0);
    if ( storage_ ) {
      value_storage_->objectiveUpdate(true);
    }
    //if ( flag ) { //&& iter>=0 ) {
      GradientSampler_->update(x);
      HessianSampler_->update(x);
      gradient_->zero();
      if ( storage_ ) {
        gradient_storage_->objectiveUpdate(true);
      }
    //}
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    initialize(x);
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

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    initialize(x);
    g.zero(); pointDual_->zero(); sumDual_->zero();
    std::vector<Ptr<Vector<Real>>> ptgs;
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

  void hessVec( Vector<Real> &hv, const Vector<Real> &v,
          const Vector<Real> &x, Real &tol ) {
    initialize(x);
    hv.zero(); pointDual_->zero(); sumDual_->zero();
    for ( int i = 0; i < HessianSampler_->numMySamples(); ++i ) {
      getHessVec(*pointDual_,v,x,HessianSampler_->getMyPoint(i),tol);
      sumDual_->axpy(HessianSampler_->getMyWeight(i),*pointDual_);
    }
    HessianSampler_->sumAll(*sumDual_,hv);
  }

  void precond( Vector<Real> &Pv, const Vector<Real> &v,
                        const Vector<Real> &x, Real &tol ) {
    Pv.set(v.dual());
  }
};

}

#endif
