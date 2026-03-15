// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_OBJECTIVEVECTOR_DEF_HPP
#define ROL_OED_OBJECTIVEVECTOR_DEF_HPP

#include "ROL_AffineTransformObjective.hpp"
#include "ROL_OED_BilinearConstraint.hpp"
#include "ROL_OED_LinearObjective.hpp"
#include "ROL_OED_QuadraticObjective.hpp"
#include "ROL_OED_A_HomObjective.hpp"
#include "ROL_OED_C_HomObjective.hpp"
#include "ROL_OED_D_HomObjective.hpp"
#include "ROL_OED_I_HomObjective.hpp"
#include "ROL_OED_Itrace_HomObjective.hpp"
#include "ROL_OED_A_HetObjective.hpp"
#include "ROL_OED_C_HetObjective.hpp"
#include "ROL_OED_D_HetObjective.hpp"
#include "ROL_OED_I_HetObjective.hpp"
#include "ROL_OED_Itrace_HetObjective.hpp"
#include "ROL_OED_Radamacher.hpp"

namespace ROL {
namespace OED {

template<typename Real>
ObjectiveArray<Real>::ObjectiveArray(const Ptr<Constraint<Real>>      &model,
                                     const Ptr<Vector<Real>>          &obs,
                                     const Ptr<SampleGenerator<Real>> &sampler,
                                     const Ptr<MomentOperator<Real>>  &cov,
                                           ParameterList              &plist)
  : modelC_(model), theta_(nullPtr), obs_(obs), cov_(cov), sampler_(sampler),
    useStorage_(plist.sublist("OED").get("Use Storage",true)),
    nfact_(0), usePF_(false) {
  RegressionType regType;
  Ptr<Noise<Real>> noise;
  cov->getRegressionInfo(regType,isHom_,noise);
}

template<typename Real>
ObjectiveArray<Real>::ObjectiveArray(const Ptr<Objective<Real>>       &model,
                                     const Ptr<SampleGenerator<Real>> &sampler,
                                     const Ptr<MomentOperator<Real>>  &cov,
                                           ParameterList              &plist)
  : modelO_(model), theta_(nullPtr), cov_(cov), sampler_(sampler),
    useStorage_(plist.sublist("OED").get("Use Storage",true)),
    nfact_(0), usePF_(false) {
  RegressionType regType;
  Ptr<Noise<Real>> noise;
  cov->getRegressionInfo(regType,isHom_,noise);
}

template<typename Real>
void ObjectiveArray<Real>::setProbabilityVector(const Ptr<Vector<Real>> &p) {
  if (theta_==nullPtr) {
    PFop_  = makePtr<DiagonalOperator<Real>>(*p);
    PFvc_  = p->clone(); PFvc_->zero();
    usePF_ = true;
  }
}

template<typename Real>
void ObjectiveArray<Real>::setTheta(const Ptr<Vector<Real>> &theta) {
  usePF_ = false;
  theta_ = theta;
}

template<typename Real>
Ptr<Objective<Real>> ObjectiveArray<Real>::buildHomObjective(ParameterList &plist,
                                                             const Ptr<Factors<Real>> &factors,
                                                             const Ptr<MomentOperator<Real>> &cov0,
                                                             const Ptr<Vector<Real>> &theta,
                                                             const Ptr<SampleGenerator<Real>> &sampler,
                                                             const Ptr<Objective<Real>> &predFun) {
  Ptr<Vector<Real>> theta0 = theta;
  if (theta_ != nullPtr) theta0 = theta_;

  type_ = plist.sublist("OED").get("Optimality Type","C");
  Ptr<Objective<Real>> obj, sobj;
  bool useTrace = plist.sublist("OED").sublist("I-Optimality").get("Use Trace Form",false);
  if ((type_ == "I" && !useTrace) || type_ == "R") {
    Ptr<BilinearConstraint<Real>> covar0;
    Ptr<LinearObjective<Real>> lobj;
    if (predFun == nullPtr) {
      covar0 = makePtr<BilinearConstraint<Real>>(factors,cov0,type_);
      lobj   = makePtr<LinearObjective<Real>>(factors,type_);
    }
    else {
      auto factorsPV = makePtr<Factors<Real>>(predFun,theta0,sampler);
      covar0 = makePtr<BilinearConstraint<Real>>(factorsPV,cov0,type_);
      lobj   = makePtr<LinearObjective<Real>>(factorsPV,type_);
    }
    obj     = makePtr<Hom::I_Objective<Real>>(covar0,lobj,theta0,useStorage_);
    int dim = sampler_->getMyPoint(0).size();
    Real wt = 0.0;
    std::vector<Real> sum(dim,0.0), mean(dim,0.0), pt(dim,0.0);
    for (int i = 0; i < sampler_->numMySamples(); ++i) {
      pt = sampler_->getMyPoint(i);
      wt = sampler_->getMyWeight(i);
      for (int j = 0; j < dim; ++j) sum[j] += wt*pt[j];
    }
    sampler_->sumAll(&sum[0],&mean[0],dim);
    obj->setParameter(mean);
  }
  else if (type_ == "C") {
    Ptr<Vector<Real>> c = theta0->dual().clone();
    Real cval = plist.sublist("OED").sublist("C-Optimality").get("C Value",1.0);
    c->setScalar(cval);
    auto covar0 = makePtr<BilinearConstraint<Real>>(factors,cov0,c);
    auto lobj   = makePtr<LinearObjective<Real>>(c);
    obj         = makePtr<Hom::C_Objective<Real>>(covar0,lobj,theta0,useStorage_);
  }
  else if (type_ == "D") {
    auto traceSampler = makePtr<TraceSampler<Real>>(theta0);
    auto covar0       = makePtr<BilinearConstraint<Real>>(factors,cov0,type_,traceSampler);
    obj               = makePtr<Hom::D_Objective<Real>>(covar0,theta0,useStorage_);
  }
  else if (type_ == "A") {
    bool useRandomTrace = plist.sublist("OED").sublist("A-Optimality").get("Randomized Trace Estimation",false);
    int nRandomTrace    = plist.sublist("OED").sublist("A-Optimality").get("Number of Samples",100);
    int nfactors        = theta0->dimension();
    int size            = (useRandomTrace ? nRandomTrace : nfactors);
    Ptr<TraceSampler<Real>> traceSampler;
    if (useRandomTrace) traceSampler = makePtr<Radamacher<Real>>(theta0,size);
    else                traceSampler = makePtr<TraceSampler<Real>>(theta0);
    const int one(1);
    Real val = (useRandomTrace ? one/static_cast<Real>(nRandomTrace) : one);
    std::vector<Real> weight(size,val);
    auto covar0 = makePtr<BilinearConstraint<Real>>(factors,cov0,type_,traceSampler);
    auto lobj   = makePtr<LinearObjective<Real>>(theta0,traceSampler);
    obj         = makePtr<Hom::A_Objective<Real>>(covar0,lobj,theta0,weight,useStorage_);
    obj->setParameter({static_cast<Real>(0)});
  }
  else if (type_ == "I" && useTrace) {
    bool useRandomTrace = plist.sublist("OED").sublist("I-Optimality").get("Randomized Trace Estimation",false);
    int nRandomTrace    = plist.sublist("OED").sublist("I-Optimality").get("Number of Samples",100);
    int nfactors        = theta0->dimension();
    int size            = (useRandomTrace ? nRandomTrace : nfactors);
    Ptr<TraceSampler<Real>> traceSampler;
    if (useRandomTrace) traceSampler = makePtr<Radamacher<Real>>(theta0,size);
    else                traceSampler = makePtr<TraceSampler<Real>>(theta0);
    const int one(1);
    Real val = (useRandomTrace ? one/static_cast<Real>(nRandomTrace) : one);
    std::vector<Real> weight(size,val);
    Ptr<BilinearConstraint<Real>> covar0;
    if (predFun == nullPtr) {
      covar0 = makePtr<BilinearConstraint<Real>>(factors,cov0,type_,traceSampler);
    }
    else {
      auto factorsPV = makePtr<Factors<Real>>(predFun,theta0,sampler);
      covar0 = makePtr<BilinearConstraint<Real>>(factorsPV,cov0,type_,traceSampler);
    }
    auto lobj = makePtr<LinearObjective<Real>>(theta0,traceSampler);
    obj       = makePtr<Hom::Itrace_Objective<Real>>(covar0,lobj,theta0,sampler,weight,useStorage_);
    obj->setParameter({static_cast<Real>(0)});
  }
  else {
    throw Exception::NotImplemented(">>> OED::Factory : Optimality type not implemented!");
  }
  auto obj0 = obj;
  if (theta_ != nullPtr) {
    PFop_ = makePtr<DiagonalOperator<Real>>(*theta);
    PFvc_ = theta->clone(); PFvc_->zero();
  }
  if (theta_ != nullPtr || usePF_)
    obj0 = makePtr<AffineTransformObjective<Real>>(obj,PFop_,PFvc_);
  return obj0;
}

template<typename Real>
Ptr<Objective<Real>> ObjectiveArray<Real>::buildHomObjective(const Ptr<Vector<Real>> &c,
                                                             const Ptr<Factors<Real>> &factors,
                                                             const Ptr<MomentOperator<Real>> &cov0,
                                                             const Ptr<Vector<Real>> &theta) {
  Ptr<Vector<Real>> theta0 = theta;
  if (theta_ != nullPtr) theta0 = theta_;

  auto covar0 = makePtr<BilinearConstraint<Real>>(factors,cov0,c);
  auto lobj   = makePtr<LinearObjective<Real>>(c);
  auto obj    = makePtr<Hom::C_Objective<Real>>(covar0,lobj,theta0,useStorage_);
  auto obj0   = obj;
  if (theta_ != nullPtr) {
    PFop_ = makePtr<DiagonalOperator<Real>>(*theta);
    PFvc_ = theta->clone(); PFvc_->zero();
  }
  if (theta_ != nullPtr || usePF_)
    obj0 = makePtr<AffineTransformObjective<Real>>(obj,PFop_,PFvc_);
  return obj0;
}

template<typename Real>
Ptr<Objective<Real>> ObjectiveArray<Real>::buildHetObjective(ParameterList &plist,
                                                             const Ptr<Factors<Real>> &factors,
                                                             const Ptr<MomentOperator<Real>> &cov0,
                                                             const Ptr<MomentOperator<Real>> &cov1,
                                                             const Ptr<Vector<Real>> &theta,
                                                             const Ptr<SampleGenerator<Real>> &sampler,
                                                             const Ptr<Objective<Real>> &predFun) {
  Ptr<Vector<Real>> theta0 = theta;
  if (theta_ != nullPtr) theta0 = theta_;

  type_ = plist.sublist("OED").get("Optimality Type","C");
  Ptr<Objective<Real>> obj;
  Ptr<BilinearConstraint<Real>> covar0, covar1;
  Ptr<QuadraticObjective<Real>> qobj;
  bool useTrace = plist.sublist("OED").sublist("I-Optimality").get("Use Trace Form",false);
  if ((type_ == "I" && !useTrace) || type_ == "R") {
    covar0 = makePtr<BilinearConstraint<Real>>(factors,cov0,type_);
    if (predFun == nullPtr) {
      covar1 = makePtr<BilinearConstraint<Real>>(factors,cov1,type_);
    }
    else {
      auto factorsPV = makePtr<Factors<Real>>(predFun,theta0,sampler);
      covar1 = makePtr<BilinearConstraint<Real>>(factorsPV,cov1,type_);
    }
    qobj     = makePtr<QuadraticObjective<Real>>(covar0);
    obj      = makePtr<Het::I_Objective<Real>>(covar1,qobj,theta0,useStorage_);
    int dim  = sampler_->getMyPoint(0).size();
    Real wt(0);
    std::vector<Real> sum(dim,0.0), mean(dim,0.0), pt(dim,0.0);
    for (int i = 0; i < sampler_->numMySamples(); ++i) {
      pt = sampler_->getMyPoint(i);
      wt = sampler_->getMyWeight(i);
      for (int j = 0; j < dim; ++j) {
        sum[j] += wt*pt[j];
      }
    }
    sampler_->sumAll(&sum[0],&mean[0],dim);
    obj->setParameter(mean);
  }
  else if (type_ == "C") {
    Ptr<Vector<Real>> c = theta0->dual().clone();
    Real cval = plist.sublist("OED").sublist("C-Optimality").get("C Value",1.0);
    c->setScalar(cval);
    covar0 = makePtr<BilinearConstraint<Real>>(factors,cov0,c);
    covar1 = makePtr<BilinearConstraint<Real>>(factors,cov1,c);
    qobj   = makePtr<QuadraticObjective<Real>>(covar0);
    obj    = makePtr<Het::C_Objective<Real>>(covar1,qobj,theta0,useStorage_);
  }
  else if (type_ == "D") {
    auto traceSampler = makePtr<TraceSampler<Real>>(theta0);
    covar0            = makePtr<BilinearConstraint<Real>>(factors,cov0,type_,traceSampler);
    covar1            = makePtr<BilinearConstraint<Real>>(factors,cov1,type_,traceSampler);
    obj               = makePtr<Het::D_Objective<Real>>(covar0,covar1,theta0,useStorage_);
  }
  else if (type_ == "A") {
    bool useRandomTrace = plist.sublist("OED").sublist("A-Optimality").get("Randomized Trace Estimation",false);
    int nRandomTrace    = plist.sublist("OED").sublist("A-Optimality").get("Number of Samples",100);
    int nfactors = theta0->dimension();
    int size = (useRandomTrace ? nRandomTrace : nfactors);
    Ptr<TraceSampler<Real>> traceSampler;
    if (useRandomTrace) traceSampler = makePtr<Radamacher<Real>>(theta0,size);
    else                traceSampler = makePtr<TraceSampler<Real>>(theta0);
    const int one(1);
    Real val = (useRandomTrace ? one/static_cast<Real>(nRandomTrace) : one);
    std::vector<Real> weight(size,val);
    covar0 = makePtr<BilinearConstraint<Real>>(factors,cov0,type_,traceSampler);
    covar1 = makePtr<BilinearConstraint<Real>>(factors,cov1,type_,traceSampler);
    qobj   = makePtr<QuadraticObjective<Real>>(covar0);
    obj    = makePtr<Het::A_Objective<Real>>(covar1,qobj,theta0,weight,useStorage_);
    obj->setParameter({static_cast<Real>(0)});
  }
  else if (type_ == "I" && useTrace) {
    bool useRandomTrace = plist.sublist("OED").sublist("I-Optimality").get("Randomized Trace Estimation",false);
    int nRandomTrace    = plist.sublist("OED").sublist("I-Optimality").get("Number of Samples",100);
    int nfactors = theta0->dimension();
    int size = (useRandomTrace ? nRandomTrace : nfactors);
    Ptr<TraceSampler<Real>> traceSampler;
    if (useRandomTrace) traceSampler = makePtr<Radamacher<Real>>(theta0,size);
    else                traceSampler = makePtr<TraceSampler<Real>>(theta0);
    const int one(1);
    Real val = (useRandomTrace ? one/static_cast<Real>(nRandomTrace) : one);
    std::vector<Real> weight(size,val);
    covar0 = makePtr<BilinearConstraint<Real>>(factors,cov0,type_,traceSampler);
    if (predFun == nullPtr) {
      covar1 = makePtr<BilinearConstraint<Real>>(factors,cov1,type_,traceSampler);
    }
    else {
      auto factorsPV = makePtr<Factors<Real>>(predFun,theta0,sampler);
      covar1 = makePtr<BilinearConstraint<Real>>(factorsPV,cov1,type_,traceSampler);
    }
    qobj = makePtr<QuadraticObjective<Real>>(covar0);
    obj  = makePtr<Het::Itrace_Objective<Real>>(covar1,qobj,theta0,sampler,weight,useStorage_);
    obj->setParameter({static_cast<Real>(0)});
  }
  else {
    throw Exception::NotImplemented(">>> OED::Factory : Optimality type not implemented!");
  }
  auto obj0 = obj;
  if (theta_ != nullPtr) {
    PFop_ = makePtr<DiagonalOperator<Real>>(*theta);
    PFvc_ = theta->clone(); PFvc_->zero();
  }
  if (theta_ != nullPtr || usePF_)
    obj0 = makePtr<AffineTransformObjective<Real>>(obj,PFop_,PFvc_);
  return obj0;
}

template<typename Real>
Ptr<Objective<Real>> ObjectiveArray<Real>::buildHetObjective(const Ptr<Vector<Real>> &c,
                                                             const Ptr<Factors<Real>> &factors,
                                                             const Ptr<MomentOperator<Real>> &cov0,
                                                             const Ptr<MomentOperator<Real>> &cov1,
                                                             const Ptr<Vector<Real>> &theta) {
  Ptr<Vector<Real>> theta0 = theta;
  if (theta_ != nullPtr) theta0 = theta_;

  auto covar0 = makePtr<BilinearConstraint<Real>>(factors,cov0,c);
  auto covar1 = makePtr<BilinearConstraint<Real>>(factors,cov1,c);
  auto qobj   = makePtr<QuadraticObjective<Real>>(covar0);
  auto obj    = makePtr<Het::C_Objective<Real>>(covar1,qobj,theta,useStorage_);
  auto obj0   = obj;
  if (theta_ != nullPtr) {
    PFop_ = makePtr<DiagonalOperator<Real>>(*theta);
    PFvc_ = theta->clone(); PFvc_->zero();
  }
  if (theta_ != nullPtr || usePF_)
    obj0 = makePtr<AffineTransformObjective<Real>>(obj,PFop_,PFvc_);
  return obj0;
}

template<typename Real>
void ObjectiveArray<Real>::addObjective(const Ptr<Vector<Real>>          &theta,
                                        ParameterList                    &plist,
                                        Real                              weight,
                                        const Ptr<SampleGenerator<Real>> &predSamp,
                                        const Ptr<Objective<Real>>       &predFunc) {
  Ptr<Vector<Real>> theta0 = theta;
  if (theta_ != nullPtr) theta0 = theta_;

  Ptr<MomentOperator<Real>> cov0, cov1;
  cov0 = cov_->clone();
  cov0->setMatrixNumber(0);
  Ptr<LinearOperator<Real>> P = cov_->getPerturbation();
  if (P != nullPtr)
    cov0->setPerturbation(P);
  if (modelO_ == nullPtr)
    cov0->generateFactors(modelC_,theta0,obs_,sampler_);
  else
    cov0->generateFactors(modelO_,theta0,sampler_);
  Ptr<Factors<Real>> factors = cov0->getFactors();
  if (!isHom_) {
    cov1 = cov0->clone();
    cov1->setFactors(factors);
    cov1->setMatrixNumber(1);
    if (P != nullPtr)
      cov1->setPerturbation(P);
  }
  objVec_.push_back(nullPtr);
  if (isHom_)
    objVec_[nfact_] = buildHomObjective(plist,factors,cov0,theta,predSamp,predFunc);
  else
    objVec_[nfact_] = buildHetObjective(plist,factors,cov0,cov1,theta,predSamp,predFunc);
  weights_.push_back(weight);
  nfact_++;
}

template<typename Real>
void ObjectiveArray<Real>::addObjective(const Ptr<Vector<Real>> &theta,
                                        const Ptr<Vector<Real>> &c,
                                        Real                     weight) {
  Ptr<Vector<Real>> theta0 = theta;
  if (theta_ != nullPtr) theta0 = theta_;

  Ptr<MomentOperator<Real>> cov0, cov1;
  cov0 = cov_->clone();
  cov0->setMatrixNumber(0);
  if (modelO_ == nullPtr)
    cov0->generateFactors(modelC_,theta0,obs_,sampler_);
  else
    cov0->generateFactors(modelO_,theta0,sampler_);
  Ptr<LinearOperator<Real>> P = cov_->getPerturbation();
  Ptr<Factors<Real>> factors = cov0->getFactors();
  if (P != nullPtr)
    cov0->setPerturbation(P);
  if (!isHom_) {
    cov1 = cov0->clone();
    cov1->setFactors(factors);
    cov1->setMatrixNumber(1);
    if (P != nullPtr)
      cov1->setPerturbation(P);
  }
  objVec_.push_back(nullPtr);
  if (isHom_)
    objVec_[nfact_] = buildHomObjective(c,factors,cov0,theta);
  else
    objVec_[nfact_] = buildHetObjective(c,factors,cov0,cov1,theta);
  weights_.push_back(weight);
  nfact_++;
}

template<typename Real>
const Ptr<Objective<Real>> ObjectiveArray<Real>::getObjective(unsigned k) const {
  if (k < nfact_) return objVec_[k];
  else
    throw Exception::NotImplemented(">>> OED::ObjectiveArray : k is out of range!");
  return nullPtr;
}

template<typename Real>
Real ObjectiveArray<Real>::getWeight(unsigned k) const {
  if (k < nfact_) return weights_[k];
  else
    throw Exception::NotImplemented(">>> OED::ObjectiveArray : k is out of range!");
  return static_cast<Real>(0);
}

template<typename Real>
std::vector<Real> ObjectiveArray<Real>::getWeights() const {
  return weights_;
}

template<typename Real>
const Ptr<Factors<Real>> ObjectiveArray<Real>::getFactors(const Ptr<Vector<Real>> &theta) const {
  auto cov = cov_->clone();
  if (modelO_ == nullPtr)
    cov->generateFactors(modelC_,theta,obs_,sampler_);
  else
    cov->generateFactors(modelO_,theta,sampler_);
  return cov->getFactors();
  //if (modelO_ == nullPtr)
  //  return makePtr<Factors<Real>>(modelC_,theta,obs_,sampler_);
  //else
  //  return makePtr<Factors<Real>>(modelO_,theta,sampler_);
}

template<typename Real>
const Ptr<MomentOperator<Real>> ObjectiveArray<Real>::getBaseMomentOperator() const {
  return cov_;
}

} // End OED Namespace
} // End ROL Namespace

#endif
