// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STOCHASTICOBJECTIVE_HPP
#define ROL_STOCHASTICOBJECTIVE_HPP

#include "ROL_RandVarFunctional.hpp"
#include "ROL_Objective.hpp"
#include "ROL_RiskVector.hpp"
#include "ROL_RandVarFunctionalFactory.hpp"
#include "ROL_ConvexCombinationRiskMeasure.hpp"

namespace ROL {

template<class Real>
class StochasticObjective : public Objective<Real> {
private:
  // Objective function definition
  Ptr<Objective<Real>>         obj_; // Uncertain objective function
  Ptr<RandVarFunctional<Real>> rvf_; // Random variable functional

  // Sampler generators
  Ptr<SampleGenerator<Real>>  vsampler_; // Sampler for objective value
  Ptr<SampleGenerator<Real>>  gsampler_; // Sampler for objective gradient
  Ptr<SampleGenerator<Real>>  hsampler_; // Sampler for objective Hessian-times-a-vector

  const int comp_;
  int index_;

  Ptr<const Vector<Real>> getConstVector(const Vector<Real> &x) const {
    const RiskVector<Real> &xrv = dynamic_cast<const RiskVector<Real>&>(x);
    return xrv.getVector();
  }

  Ptr<Vector<Real>> getVector(Vector<Real> &x) const {
    RiskVector<Real> &xrv = dynamic_cast<RiskVector<Real>&>(x);
    return xrv.getVector();
  }

  Ptr<const std::vector<Real>> getConstStat(const Vector<Real> &x) const {
    const RiskVector<Real> &xrv = dynamic_cast<const RiskVector<Real>&>(x);
    Ptr<const std::vector<Real>> xstat = xrv.getStatistic(comp_,index_);
    if (xstat == nullPtr) {
      xstat = makePtr<const std::vector<Real>>(0);
    }
    return xstat;
  }

  Ptr<std::vector<Real>> getStat(Vector<Real> &x) const {
    RiskVector<Real> &xrv = dynamic_cast<RiskVector<Real>&>(x);
    Ptr<std::vector<Real>> xstat = xrv.getStatistic(comp_,index_);
    if (xstat == nullPtr) {
      xstat = makePtr<std::vector<Real>>(0);
    }
    return xstat;
  }

public:
  virtual ~StochasticObjective() {}

  StochasticObjective( const Ptr<Objective<Real>> &obj,
               const Ptr<RandVarFunctional<Real>> &rvf,
               const Ptr<SampleGenerator<Real>>   &vsampler,
               const Ptr<SampleGenerator<Real>>   &gsampler,
               const Ptr<SampleGenerator<Real>>   &hsampler,
               const bool storage = true,
               const int comp = 0, const int index = 0 )
    : obj_(obj), rvf_(rvf),
      vsampler_(vsampler), gsampler_(gsampler), hsampler_(hsampler),
      comp_(comp), index_(index) {
    rvf->useStorage(storage);
  }

  StochasticObjective( const Ptr<Objective<Real>> &obj,
               const Ptr<RandVarFunctional<Real>> &rvf,
               const Ptr<SampleGenerator<Real>>   &vsampler,
               const Ptr<SampleGenerator<Real>>   &gsampler,
               const bool storage = true,
               const int comp = 0, const int index = 0 )
    : StochasticObjective(obj,rvf,vsampler,gsampler,gsampler,storage,comp,index) {}

  StochasticObjective( const Ptr<Objective<Real>> &obj,
               const Ptr<RandVarFunctional<Real>> &rvf,
               const Ptr<SampleGenerator<Real>>   &sampler,
               const bool storage = true,
               const int comp = 0, const int index = 0 )
    : StochasticObjective(obj,rvf,sampler,sampler,sampler,storage,comp,index) {}

  StochasticObjective( const Ptr<Objective<Real>> &obj,
               ParameterList                      &parlist,
               const Ptr<SampleGenerator<Real>>   &vsampler,
               const Ptr<SampleGenerator<Real>>   &gsampler,
               const Ptr<SampleGenerator<Real>>   &hsampler,
               const int comp = 0, const int index = 0 )
    : obj_(obj),
      vsampler_(vsampler), gsampler_(gsampler), hsampler_(hsampler),
      comp_(comp), index_(index) {
    std::string name, type = parlist.sublist("SOL").get("Type","Risk Averse");
    if (type == "Risk Averse")
      name = parlist.sublist("SOL").sublist("Risk Measure").get("Name","CVaR");

    if (type == "Risk Averse" && name == "Convex Combination Risk Measure")
      rvf_ = makePtr<ConvexCombinationRiskMeasure<Real>>(parlist);
    else
      rvf_ = RandVarFunctionalFactory<Real>(parlist);

    bool storage = parlist.sublist("SOL").get("Store Sampled Value and Gradient",true);
    rvf_->useStorage(storage);
  }

  StochasticObjective( const Ptr<Objective<Real>> &obj,
               ROL::ParameterList             &parlist,
               const Ptr<SampleGenerator<Real>>   &vsampler,
               const Ptr<SampleGenerator<Real>>   &gsampler,
               const int comp = 0, const int index = 0 )
    : StochasticObjective(obj,parlist,vsampler,gsampler,gsampler,comp,index) {}

  StochasticObjective( const Ptr<Objective<Real>> &obj,
               ROL::ParameterList             &parlist,
               const Ptr<SampleGenerator<Real>>   &sampler, 
               const int comp = 0, const int index = 0 )
    : StochasticObjective(obj,parlist,sampler,sampler,sampler,comp,index) {}

  Real computeStatistic(const Vector<Real> &x) const {
    Ptr<const std::vector<Real>> xstat = getConstStat(x);
    return rvf_->computeStatistic(xstat);
  }

  void setIndex(int ind) {
    index_ = ind;
  }

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) {
    Ptr<const Vector<Real>> x0 = getConstVector(x);
    // Update random variable functional
    rvf_->resetStorage(type);
    // Update uncertain objective function
    obj_->update(*x0,type,iter);
    // Update samplers
    vsampler_->update(*x0);
    if ( type != UpdateType::Trial && type != UpdateType::Revert ) {
      gsampler_->update(*x0);
      hsampler_->update(*x0);
    }
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    Ptr<const Vector<Real>> x0 = getConstVector(x);
    // Update random variable functional
    rvf_->resetStorage(flag);
    // Update uncertain objective function
    obj_->update(*x0,flag,iter);
    // Update samplers
    vsampler_->update(*x0);
    if ( flag ) {
      gsampler_->update(*x0);
      hsampler_->update(*x0);
    }
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    Ptr<const Vector<Real>>      x0    = getConstVector(x);
    Ptr<const std::vector<Real>> xstat = getConstStat(x);
    rvf_->initialize(*x0);
    Real val(0);
    for ( int i = 0; i < vsampler_->numMySamples(); i++ ) {
      rvf_->setSample(vsampler_->getMyPoint(i),vsampler_->getMyWeight(i));
      rvf_->updateValue(*obj_,*x0,*xstat,tol);
    }
    val = rvf_->getValue(*x0,*xstat,*vsampler_);
    return val;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    g.zero();
    Ptr<const Vector<Real>>      x0    = getConstVector(x);
    Ptr<const std::vector<Real>> xstat = getConstStat(x);
    Ptr<Vector<Real>>            g0    = getVector(g);
    Ptr<std::vector<Real>>       gstat = getStat(g);
    rvf_->initialize(*x0);
    for ( int i = 0; i < gsampler_->numMySamples(); i++ ) {
      rvf_->setSample(gsampler_->getMyPoint(i),gsampler_->getMyWeight(i));
      rvf_->updateGradient(*obj_,*x0,*xstat,tol);
    }
    rvf_->getGradient(*g0,*gstat,*x0,*xstat,*gsampler_);
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v,
                const Vector<Real> &x, Real &tol ) {
    hv.zero();
    Ptr<const Vector<Real>>      x0     = getConstVector(x);
    Ptr<const std::vector<Real>> xstat  = getConstStat(x);
    Ptr<const Vector<Real>>      v0     = getConstVector(v);
    Ptr<const std::vector<Real>> vstat  = getConstStat(v);
    Ptr<Vector<Real>>            hv0    = getVector(hv);
    Ptr<std::vector<Real>>       hvstat = getStat(hv);
    rvf_->initialize(*x0);
    for ( int i = 0; i < hsampler_->numMySamples(); i++ ) {
      rvf_->setSample(hsampler_->getMyPoint(i),hsampler_->getMyWeight(i));
      rvf_->updateHessVec(*obj_,*v0,*vstat,*x0,*xstat,tol);
    }
    rvf_->getHessVec(*hv0,*hvstat,*v0,*vstat,*x0,*xstat,*hsampler_);
  }

  virtual void precond( Vector<Real> &Pv, const Vector<Real> &v,
                  const Vector<Real> &x, Real &tol ) {
    Pv.set(v.dual());
  }
};

}

#endif
