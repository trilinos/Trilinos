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

namespace ROL {
namespace OED {

template<typename Real>
IndependentFactors<Real>::IndependentFactors(const Ptr<Objective<Real>>       &model,
                                             const Ptr<Vector<Real>>          &theta,
                                             const Ptr<SampleGenerator<Real>> &sampler,
                                             bool                              storage,
                                             bool                              ortho) {
  auto pmodel = staticPtrCast<SeparableObjective<Real>>(model);
  auto pvec   = staticPtrCast<PartitionedVector<Real>>(theta);
  numOp_ = pvec->numVectors();
  factors_.resize(numOp_);
  for (unsigned i=0; i < numOp_; ++i)
    factors_[i] = makePtr<Factors<Real>>(pmodel.get(i),pvec.get(i),sampler,storage,ortho);
}

template<typename Real>
IndependentFactors<Real>::IndependentFactors(const std::vector<Ptr<Factors<Real>>> &factors)
  : factors_(factors), numOp_(factors.size()) {
  Fk_.clear(); Fk_.resize(numOp_);
  for (unsigned i = 0u; i < numOp_; ++i)
    Fk_[i] = factors[i]->get(0)->clone();
}

template<typename Real>
void IndependentFactors<Real>::setPredictionVector(const Vector<Real> &c) {
  Factors<Real>::setPredictionVector(c);
  for (unsigned i = 0u; i < numOp_; ++i) {
    Ptr<const Vector<Real>> ci = dynamic_cast<const PartitionedVector<Real>&>(c).get(i);
    factors_[i]->setPredictionVector(*ci);
  }
}

template<typename Real>
void IndependentFactors<Real>::getPredictionVector(Vector<Real> &c) const {
  Factors<Real>::getPredictionVector(c);
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

// Compute c^T F[k] x
template<typename Real>
Real IndependentFactors<Real>::apply(const Vector<Real> &x, int k) const {
  //startTimer("apply");
  Real val(0);
  for (unsigned i = 0u; i < numOp_; ++i)
    val += factors_[i]->apply(x,k);
  //stopTimer("apply");
  return val;
}

// Compute Fx = F[k] x
template<typename Real>
void IndependentFactors<Real>::apply(Vector<Real> &Fx, const Vector<Real> &x, int k) const {
  //startTimer("apply");
  for (unsigned i = 0u; i < numOp_; ++i) {
    Ptr<Vector<Real>> Fxi = dynamic_cast<PartitionedVector<Real>&>(Fx).get(i);
    Ptr<const Vector<Real>> xi = dynamic_cast<const PartitionedVector<Real>&>(x).get(i);
    factors_[i]->apply(*Fxi,*xi,k);
  }
  //stopTimer("apply");
}

// Compute Mx = F[k]^T R F[k] x
template<typename Real>
void IndependentFactors<Real>::applyProduct(Vector<Real> &Mx, const Vector<Real> &x, int k) const {
  //startTimer("applyProduct");
  for (unsigned i = 0u; i < numOp_; ++i) {
    Ptr<Vector<Real>> Mxi = dynamic_cast<PartitionedVector<Real>&>(Mx).get(i);
    Ptr<const Vector<Real>> xi = dynamic_cast<const PartitionedVector<Real>&>(x).get(i);
    factors_[i]->applyProduct(*Mxi,*xi,k);
  }
  //stopTimer("applyProduct");
}

// Compute y^T F[k]^T R F[k] x
template<typename Real>
Real IndependentFactors<Real>::applyProduct2(const Vector<Real> &x, const Vector<Real> &y, int k) const {
  //startTimer("applyProduct2");
  Real val(0);
  for (unsigned i = 0u; i < numOp_; ++i)
    val += factors_[i]->applyProduct2(x,y,k);
  //stopTimer("applyProduct2");
  return val;
}

// Get F[k]^T c
template<typename Real>
const Ptr<const Vector<Real>> IndependentFactors<Real>::get(int k) const {
  //startTimer("get");
  for (unsigned i = 0u; i < numOp_; ++i)
    Fk_[i]->set(*factors_[i]->get(k));
  //stopTimer("get");
  return makePtr<PartitionedVector<Real>>(Fk_);
}

// Compute F(param)^T c
template<typename Real>
void IndependentFactors<Real>::evaluate(Vector<Real> &F, const std::vector<Real> &param) const {
  //startTimer("evaluate");
  for (unsigned i = 0u; i < numOp_; ++i) {
    Ptr<Vector<Real>> Fi = dynamic_cast<PartitionedVector<Real>&>(F).get(i);
    factors_[i]->evaluate(*Fi,param);
  }
  //stopTimer("evaluate");
}

template<typename Real>
int IndependentFactors<Real>::numFactors() const {
  return factors_[0]->numFactors();
}

template<typename Real>
int IndependentFactors<Real>::numMySamples() const {
  return factors_[0]->numMySamples();
}

template<typename Real>
std::vector<Real> IndependentFactors<Real>::getSample(int k) const {
  return factors_[0]->getSample(k);
}

template<typename Real>
const Ptr<Factors<Real>> IndependentFactors<Real>::getFactors(unsigned i) const {
	std::cout << i << std::endl;
  return factors_[i];
}

} // End OED Namespace
} // End ROL Namespace

#endif
