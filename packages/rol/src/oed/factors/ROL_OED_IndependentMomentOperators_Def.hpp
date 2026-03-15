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

namespace ROL {
namespace OED {

template<typename Real>
IndependentMomentOperator<Real>::IndependentMomentOperator(const std::vector<const Ptr<MomentOperator<Real>>> & Mvec)
  : ProfiledClass<Real,std::string>("OED::IndependentMomentOperator"), Mvec_(Mvec) {}

template<typename Real>
Ptr<MomentOperator<Real>> IndependentMomentOperator<Real>::clone() const {
  return makePtr<IndependentMomentOperator<Real>>(Mvec);
}

template<typename Real>
void IndependentMomentOperator<Real>::setMatrixNumber(int matNum) {
  int size = Mvec.size();
  for (int i = 0; i < size; ++i)
    Mvec_->setMatrixNumber(matNum);
}

template<typename Real>
void IndependentMomentOperator<Real>::setFactors(const Ptr<Factors<Real>> &factors) {}

template<typename Real>
void IndependentMomentOperator<Real>::setPerturbation(const Ptr<LinearOperator<Real>> &pOp) {}

template<typename Real>
void IndependentMomentOperator<Real>::apply(Vector<Real> &Mx,
                                      const Vector<Real> &x,
                                      const Vector<Real> &p) {
  
  startTimer("apply");
  int size = Mvec.size();
  for (int i = 0; i < size; ++i) {
    Ptr<Vector<Real>>       Mxi = dynamic_cast<PartitionedVector<Real>&>(Mx).get(i);
    Ptr<const Vector<Real>>  xi = dynamic_cast<const PartitionedVector<Real>&>(x).get(i);
    Mvec_[i]->apply(Mxi,xi,p);
  }
  stopTimer("apply");
}

template<typename Real>
void IndependentMomentOperator<Real>::applyDeriv(Vector<Real> &Mx,
                                           const Vector<Real> &x,
                                           const Vector<Real> &q) {
  startTimer("applyDeriv");
  int size = Mvec.size();
  for (int i = 0; i < size; ++i) {
    Ptr<Vector<Real>>       Mxi = dynamic_cast<PartitionedVector<Real>&>(Mx).get(i);
    Ptr<const Vector<Real>>  xi = dynamic_cast<const PartitionedVector<Real>&>(x).get(i);
    Mvec_[i]->applyDeriv(Mxi,xi,q);
  }
  stopTimer("applyDeriv");
}

template<typename Real>
void IndependentMomentOperator<Real>::applyInverse(Vector<Real> &Mx,
                                             const Vector<Real> &x,
                                             const Vector<Real> &p) {
  startTimer("applyInverse");
  int size = Mvec.size();
  for (int i = 0; i < size; ++i) {
    Ptr<Vector<Real>>       Mxi = dynamic_cast<PartitionedVector<Real>&>(Mx).get(i);
    Ptr<const Vector<Real>>  xi = dynamic_cast<const PartitionedVector<Real>&>(x).get(i);
    Mvec_[i]->applyInverse(Mxi,xi,p);
  }
  stopTimer("applyInverse");
}

template<typename Real>
void IndependentMomentOperator<Real>::applySampleMatrices(Vector<Real> &uXv, const Vector<Real> &u, const Vector<Real> &v) {
  startTimer("applySampleMatrices");
  int size = Mvec.size();
  for (int i = 0; i < size; ++i) {
    Ptr<Vector<Real>>       UXv = dynamic_cast<PartitionedVector<Real>&>(uXv).get(i);
    Ptr<const Vector<Real>>  ui = dynamic_cast<const PartitionedVector<Real>&>(u).get(i);
    Ptr<const Vector<Real>>  vi = dynamic_cast<const PartitionedVector<Real>&>(v).get(i);
    Mvec_[i]->applyInverse(Mxi,xi,p);
  }
  stopTimer("applySampleMatrices");
}

template<typename Real>
void IndependentMomentOperator<Real>::setNoise(const Ptr<Noise<Real>> &noise, bool isHom) {}

template<typename Real>
Real IndependentMomentOperator<Real>::getNoise(int k) const {
  Real val(1);
  return val;
}

template<typename Real>
void IndependentMomentOperator<Real>::getRegressionInfo(RegressionType &regType, bool &homNoise,
                       Ptr<Noise<Real>> &noise) const {
  Mvec_[0]->getRegressionInfo(regType,homNoise,noise);
}

template<typename Real>
Real IndependentMomentOperator<Real>::logDeterminant(const Vector<Real> &z) {
  startTimer("logDeterminant");
  int size = Mvec.size();
  Real val(0);
  for (int i = 0; i < size; ++i)
    val += Mvec_->logDeterminant(z);
  stopTimer("logDeterminant");
  return val;
}

} // End OED Namespace
} // End ROL Namespace

#endif
