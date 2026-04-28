// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_INDEPENDENT_FACTORS_DEF_HPP
#define ROL_OED_INDEPENDENT_FACTORS_DEF_HPP

#include "ROL_PartitionedVector.hpp"

namespace ROL::OED {

template<typename Real>
IndependentFactors<Real>::IndependentFactors(const Ptr<Objective<Real>>& model,
                                             const Ptr<const Vector<Real>>& theta,
                                             const Ptr<SampleGenerator<Real>>& sampler) {
  auto pmodel = staticPtrCast<SeparableObjective<Real>>(model);
  auto pvec   = staticPtrCast<PartitionedVector<Real>>(theta);
  numOp_ = pvec->numVectors();
  factors_.resize(numOp_);
  nobs_ = 0u;
  nfact_ = 0u;
  for (unsigned i=0; i < numOp_; ++i) {
    factors_[i] = makePtr<Factors<Real>>(pmodel.get(i),pvec.get(i),sampler);
    nobs_ += factors_[i]->numObservations();
    nfact_ += factors_[i]->numFactors();
  }
}

template<typename Real>
IndependentFactors<Real>::IndependentFactors(const std::vector<Ptr<Factors<Real>>>& factors)
  : factors_(factors), numOp_(factors.size()) {
  nobs_ = 0u;
  nfact_ = 0u;
  for (unsigned i = 0u; i < numOp_; ++i) {
    nobs_ += factors_[i]->numObservations();
    nfact_ += factors_[i]->numFactors();
  }
}

// Create a vector in the parameter space
template<typename Real>
Ptr<Vector<Real>> IndependentFactors<Real>::createParameterVector(bool dual) const {
  std::vector<Ptr<Vector<Real>>> vecs(numOp_);
  for (unsigned i = 0; i < numOp_; ++i)
    vecs[i] = factors_[i]->createParameterVector(dual);
  return makePtr<PartitionedVector<Real>>(vecs);
}

// Create a vector in the observation space
template<typename Real>
Ptr<Vector<Real>> IndependentFactors<Real>::createObservationVector(bool dual) const {
  std::vector<Ptr<Vector<Real>>> vecs(numOp_);
  for (unsigned i = 0; i < numOp_; ++i)
    vecs[i] = factors_[i]->createObservationVector(dual);
  return makePtr<PartitionedVector<Real>>(vecs);
}

// Compute Fx = F[k] x
template<typename Real>
void IndependentFactors<Real>::apply(Vector<Real> &Fx, const Vector<Real> &x, int k) const {
  apply(Fx,x,getSample(k));
}

// Compute Fx = F[k] x
template<typename Real>
void IndependentFactors<Real>::apply(Vector<Real> &Fx, const Vector<Real> &x, const std::vector<Real>& pt) const {
  //startTimer("apply");
  for (unsigned i = 0u; i < numOp_; ++i) {
    auto Fxi = dynamic_cast<PartitionedVector<Real>&>(Fx).get(i);
    auto xi = dynamic_cast<const PartitionedVector<Real>&>(x).get(i);
    factors_[i]->apply(*Fxi,*xi,pt);
  }
  //stopTimer("apply");
}

// Compute Fx = F[k]* x
template<typename Real>
void IndependentFactors<Real>::applyAdjoint(Vector<Real> &Fx, const Vector<Real> &x, int k) const {
  apply(Fx,x,getSample(k));
}

// Compute Fx = F[k]* x
template<typename Real>
void IndependentFactors<Real>::applyAdjoint(Vector<Real> &Fx, const Vector<Real> &x, const std::vector<Real>& pt) const {
  //startTimer("applyAdjoint");
  for (unsigned i = 0u; i < numOp_; ++i) {
    auto Fxi = dynamic_cast<PartitionedVector<Real>&>(Fx).get(i);
    auto xi = dynamic_cast<const PartitionedVector<Real>&>(x).get(i);
    factors_[i]->applyAdjoint(*Fxi,*xi,pt);
  }
  //stopTimer("applyAdjoint");
}

} // End ROL::OED Namespace

#endif
