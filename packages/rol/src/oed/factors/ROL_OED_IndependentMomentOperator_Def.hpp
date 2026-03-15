// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_INDEPENDENT_COVARIANCE_OPERATOR_DEF_HPP
#define ROL_OED_INDEPENDENT_COVARIANCE_OPERATOR_DEF_HPP

#include "ROL_PartitionedVector.hpp"
#include "ROL_BlockDiagonalOperator.hpp"
#include "ROL_SeparableObjective.hpp"
#include "ROL_SeparableConstraint.hpp"
#include "ROL_OED_IndependentFactors.hpp"

namespace ROL {
namespace OED {

template<typename Real>
IndependentMomentOperator<Real>::IndependentMomentOperator(const std::vector<Ptr<MomentOperator<Real>>> &Mvec)
  : Mvec_(Mvec), numOp_(Mvec.size()) {}
  //: ProfiledClass<Real,std::string>("OED::IndependentMomentOperator"),

template<typename Real>
void IndependentMomentOperator<Real>::update(const Vector<Real> &p, UpdateType type, int iter) {
  for (unsigned i = 0u; i < numOp_; ++i) Mvec_[i]->update(p,type,iter);
}

template<typename Real>
Ptr<MomentOperator<Real>> IndependentMomentOperator<Real>::clone() const {
  std::vector<Ptr<MomentOperator<Real>>> Mvec(numOp_);
  for (unsigned i = 0u; i < numOp_; ++i) Mvec[i] = Mvec_[i]->clone();
  return makePtr<IndependentMomentOperator<Real>>(Mvec);
}

template<typename Real>
void IndependentMomentOperator<Real>::setMatrixNumber(int matNum) {
  MomentOperator<Real>::setMatrixNumber(matNum);
  for (unsigned i = 0u; i < numOp_; ++i)
    Mvec_[i]->setMatrixNumber(matNum);
}

template<typename Real>
void IndependentMomentOperator<Real>::setFactors(const Ptr<Factors<Real>> &factors) {
  MomentOperator<Real>::setFactors(factors);
  auto ifactors = staticPtrCast<IndependentFactors<Real>>(factors);
  for (unsigned i = 0u; i < numOp_; ++i)
    Mvec_[i]->setFactors(ifactors->getFactors(i));
}

template<typename Real>
void IndependentMomentOperator<Real>::generateFactors(const Ptr<Constraint<Real>>      &model,
                                                      const Ptr<Vector<Real>>          &theta,
                                                      const Ptr<Vector<Real>>          &obs,
                                                      const Ptr<SampleGenerator<Real>> &sampler,
                                                      bool                              storage,
                                                      const Ptr<Vector<Real>>          &c,
                                                      bool                              ortho) {
  auto pmodel = staticPtrCast<SeparableConstraint<Real>>(model);
  auto ptheta = staticPtrCast<PartitionedVector<Real>>(theta);
  auto pobs   = staticPtrCast<PartitionedVector<Real>>(obs);
  std::vector<Ptr<Factors<Real>>> factors(numOp_);
  for (unsigned i = 0u; i < numOp_; ++i) {
    Ptr<Vector<Real>> ci = nullPtr;
    if (c!=nullPtr) ci = staticPtrCast<PartitionedVector<Real>>(c)->get(i);
    Mvec_[i]->generateFactors(pmodel->get(i),ptheta->get(i),pobs->get(i),sampler,storage,ci,ortho);
    factors[i] = Mvec_[i]->getFactors();
  }
  MomentOperator<Real>::setFactors(makePtr<IndependentFactors<Real>>(factors));
}

template<typename Real>
void IndependentMomentOperator<Real>::generateFactors(const Ptr<Objective<Real>>       &model,
                                                      const Ptr<Vector<Real>>          &theta,
                                                      const Ptr<SampleGenerator<Real>> &sampler,
                                                      bool                              storage,
                                                      bool                              ortho) {
  auto pmodel = staticPtrCast<SeparableObjective<Real>>(model);
  auto ptheta = staticPtrCast<PartitionedVector<Real>>(theta);
  std::vector<Ptr<Factors<Real>>> factors(numOp_);
  for (unsigned i = 0u; i < numOp_; ++i) {
    Mvec_[i]->generateFactors(pmodel->get(i),ptheta->get(i),sampler,storage,ortho);
    factors[i] = Mvec_[i]->getFactors();
  }
  MomentOperator<Real>::setFactors(makePtr<IndependentFactors<Real>>(factors));
}

template<typename Real>
void IndependentMomentOperator<Real>::setPerturbation(const Ptr<LinearOperator<Real>> &pOp) {
  MomentOperator<Real>::setPerturbation(pOp);
  auto iop = staticPtrCast<BlockDiagonalOperator<Real>>(pOp);
  for (unsigned i = 0u; i < numOp_; ++i)
    Mvec_[i]->setPerturbation(iop->get(i));
}

template<typename Real>
void IndependentMomentOperator<Real>::apply(Vector<Real> &Mx,
                                      const Vector<Real> &x,
                                      const Vector<Real> &p) {
  startTimer("apply");
  for (unsigned i = 0u; i < numOp_; ++i) {
    Ptr<Vector<Real>>       Mxi = dynamic_cast<PartitionedVector<Real>&>(Mx).get(i);
    Ptr<const Vector<Real>>  xi = dynamic_cast<const PartitionedVector<Real>&>(x).get(i);
    Mvec_[i]->apply(*Mxi,*xi,p);
  }
  stopTimer("apply");
}

template<typename Real>
void IndependentMomentOperator<Real>::applyDeriv(Vector<Real> &Mx,
                                           const Vector<Real> &x,
                                           const Vector<Real> &q) {
  startTimer("applyDeriv");
  for (unsigned i = 0u; i < numOp_; ++i) {
    Ptr<Vector<Real>>       Mxi = dynamic_cast<PartitionedVector<Real>&>(Mx).get(i);
    Ptr<const Vector<Real>>  xi = dynamic_cast<const PartitionedVector<Real>&>(x).get(i);
    Mvec_[i]->applyDeriv(*Mxi,*xi,q);
  }
  stopTimer("applyDeriv");
}

template<typename Real>
void IndependentMomentOperator<Real>::applyInverse(Vector<Real> &Mx,
                                             const Vector<Real> &x,
                                             const Vector<Real> &p) {
  startTimer("applyInverse");
  for (unsigned i = 0u; i < numOp_; ++i) {
    Ptr<Vector<Real>>       Mxi = dynamic_cast<PartitionedVector<Real>&>(Mx).get(i);
    Ptr<const Vector<Real>>  xi = dynamic_cast<const PartitionedVector<Real>&>(x).get(i);
    Mvec_[i]->applyInverse(*Mxi,*xi,p);
  }
  stopTimer("applyInverse");
}

template<typename Real>
void IndependentMomentOperator<Real>::applySampleMatrices(Vector<Real> &uXv, const Vector<Real> &u, const Vector<Real> &v) {
  startTimer("applySampleMatrices");
  auto uXvi = uXv.clone();
  uXv.zero();
  for (unsigned i = 0u; i < numOp_; ++i) {
    //Ptr<Vector<Real>>      uXvi = dynamic_cast<PartitionedVector<Real>&>(uXv).get(i);
    Ptr<const Vector<Real>>  ui = dynamic_cast<const PartitionedVector<Real>&>(u).get(i);
    Ptr<const Vector<Real>>  vi = dynamic_cast<const PartitionedVector<Real>&>(v).get(i);
    Mvec_[i]->applySampleMatrices(*uXvi,*ui,*vi);
    uXv.plus(*uXvi);
  }
  stopTimer("applySampleMatrices");
}

template<typename Real>
void IndependentMomentOperator<Real>::setNoise(const Ptr<Noise<Real>> &noise, bool isHom) {
  MomentOperator<Real>::setNoise(noise,isHom);
  for (unsigned i = 0u; i < numOp_; ++i) Mvec_[i]->setNoise(noise,isHom);
}

template<typename Real>
Real IndependentMomentOperator<Real>::getNoise(int k) const {
  return Mvec_[0]->getNoise(k);
}

template<typename Real>
void IndependentMomentOperator<Real>::getRegressionInfo(RegressionType &regType, bool &homNoise,
                       Ptr<Noise<Real>> &noise) const {
  Mvec_[0]->getRegressionInfo(regType,homNoise,noise);
}

template<typename Real>
Real IndependentMomentOperator<Real>::logDeterminant(const Vector<Real> &z) {
  startTimer("logDeterminant");
  Real val(0);
  for (unsigned i = 0u; i < numOp_; ++i)
    val += Mvec_[i]->logDeterminant(z);
  stopTimer("logDeterminant");
  return val;
}

} // End OED Namespace
} // End ROL Namespace

#endif
