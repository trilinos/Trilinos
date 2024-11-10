// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_ELASTICOBJECTIVEDEF_H
#define ROL_ELASTICOBJECTIVEDEF_H


namespace ROL {

template<typename Real>
ElasticObjective<Real>::ElasticObjective(const Ptr<Objective<Real>> &obj,
                                         const Ptr<Constraint<Real>> &con,
                                         const Real penaltyParameter,
                                         const Real sigma,
                                         const Vector<Real> &dualOptVec,
                                         const Vector<Real> &primConVec,
                                         const Vector<Real> &dualConVec,
                                         ParameterList &parlist)
  : sigma_(sigma), cscale_(1) {
  alobj_ = makePtr<AugmentedLagrangianObjective<Real>>(obj,con,penaltyParameter,dualOptVec,primConVec,dualConVec,parlist);
  e_ = primConVec.clone(); e_->setScalar(static_cast<Real>(1));
  tmp_ = primConVec.clone();
}

template<typename Real>
ElasticObjective<Real>::ElasticObjective(const Ptr<Objective<Real>> &obj,
                                         const Ptr<Constraint<Real>> &con,
                                         const Real penaltyParameter,
                                         const Real sigma,
                                         const Vector<Real> &dualOptVec,
                                         const Vector<Real> &primConVec,
                                         const Vector<Real> &dualConVec,
                                         const bool scaleLagrangian,
                                         const int HessianApprox)
  : sigma_(sigma), cscale_(1) {
  alobj_ = makePtr<AugmentedLagrangianObjective<Real>>(obj,con,penaltyParameter,dualOptVec,primConVec,dualConVec,scaleLagrangian,HessianApprox);
  e_ = primConVec.clone(); e_->setScalar(static_cast<Real>(1));
  tmp_ = primConVec.clone();
}

template<typename Real>
void ElasticObjective<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {
  Ptr<const Vector<Real>> xs = dynamic_cast<const PartitionedVector<Real>&>(x).get(0);
  alobj_->update(*xs,type,iter);
}

template<typename Real>
void ElasticObjective<Real>::setScaling(const Real fscale, const Real cscale) {
  cscale_ = cscale;
  alobj_->setScaling(fscale,cscale);
}

template<typename Real>
Real ElasticObjective<Real>::value( const Vector<Real> &x, Real &tol ) {
  Ptr<const Vector<Real>> xs = dynamic_cast<const PartitionedVector<Real>&>(x).get(0);
  Ptr<const Vector<Real>> xu = dynamic_cast<const PartitionedVector<Real>&>(x).get(1);
  Ptr<const Vector<Real>> xv = dynamic_cast<const PartitionedVector<Real>&>(x).get(2);
  Real val = alobj_->value(*xs,tol);
  tmp_->set(*xu); tmp_->plus(*xv);
  Real pen = sigma_ * cscale_ * e_->dot(*tmp_);
  return val + pen;
}

template<typename Real>
void ElasticObjective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  Ptr<Vector<Real>> gs = dynamic_cast<PartitionedVector<Real>&>(g).get(0);
  Ptr<Vector<Real>> gu = dynamic_cast<PartitionedVector<Real>&>(g).get(1);
  Ptr<Vector<Real>> gv = dynamic_cast<PartitionedVector<Real>&>(g).get(2);
  Ptr<const Vector<Real>> xs = dynamic_cast<const PartitionedVector<Real>&>(x).get(0);
  alobj_->gradient(*gs,*xs,tol);
  gu->set(*e_); gu->scale(sigma_ * cscale_);
  gv->set(*e_); gv->scale(sigma_ * cscale_);
}

template<typename Real>
void ElasticObjective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  Ptr<Vector<Real>> hvs = dynamic_cast<PartitionedVector<Real>&>(hv).get(0);
  Ptr<Vector<Real>> hvu = dynamic_cast<PartitionedVector<Real>&>(hv).get(1);
  Ptr<Vector<Real>> hvv = dynamic_cast<PartitionedVector<Real>&>(hv).get(2);
  Ptr<const Vector<Real>> vs = dynamic_cast<const PartitionedVector<Real>&>(v).get(0);
  Ptr<const Vector<Real>> xs = dynamic_cast<const PartitionedVector<Real>&>(x).get(0);
  alobj_->hessVec(*hvs,*vs,*xs,tol);
  hvu->zero();
  hvv->zero();
}

template<typename Real>
Real ElasticObjective<Real>::getObjectiveValue(const Vector<Real> &x, Real &tol) {
  return alobj_->getObjectiveValue(x,tol);
}

template<typename Real>
const Ptr<const Vector<Real>> ElasticObjective<Real>::getObjectiveGradient(const Vector<Real> &x, Real &tol) {
  return alobj_->getObjectiveGradient(x,tol);
}

template<typename Real>
const Ptr<const Vector<Real>> ElasticObjective<Real>::getConstraintVec(const Vector<Real> &x, Real &tol) {
  return alobj_->getConstraintVec(x,tol);
}

template<typename Real>
int ElasticObjective<Real>::getNumberConstraintEvaluations(void) const {
  return alobj_->getNumberConstraintEvaluations();
}

template<typename Real>
int ElasticObjective<Real>::getNumberFunctionEvaluations(void) const {
  return alobj_->getNumberFunctionEvaluations();
}

template<typename Real>
int ElasticObjective<Real>::getNumberGradientEvaluations(void) const {
  return alobj_->getNumberGradientEvaluations();
}

template<typename Real>
void ElasticObjective<Real>::reset(const Vector<Real> &multiplier, Real penaltyParameter, Real sigma) {
  sigma_ = sigma;
  alobj_->reset(multiplier,penaltyParameter);
}

template<typename Real>
const Ptr<AugmentedLagrangianObjective<Real>> ElasticObjective<Real>::getAugmentedLagrangian(void) const {
  return alobj_;
}

} // namespace ROL

#endif
