// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_AFFINE_TRANSFORM_OBJECTIVE_DEF_H
#define ROL_AFFINE_TRANSFORM_OBJECTIVE_DEF_H

namespace ROL {

template<typename Real>
AffineTransformObjective<Real>::AffineTransformObjective(const Ptr<Objective<Real>>        &obj,
                                                         const Ptr<Constraint<Real>>       &acon,
                                                         const Vector<Real>                &range,
                                                         const Ptr<VectorController<Real>> &storage)
  : obj_(obj), acon_(acon), storage_(storage) {
  primal_ = range.clone();
  Av_     = range.clone();
  dual_   = range.dual().clone();
  if (storage == nullPtr) storage_ = makePtr<VectorController<Real>>();
}

template<typename Real>
AffineTransformObjective<Real>::AffineTransformObjective(const Ptr<Objective<Real>>        &obj,
                                                         const Ptr<LinearConstraint<Real>> &acon,
                                                         const Ptr<VectorController<Real>> &storage)
  : obj_(obj), acon_(acon), storage_(storage) {
  primal_ = acon->createRangeSpaceVector();
  Av_     = acon->createRangeSpaceVector();
  dual_   = primal_->dual().clone();
  if (storage == nullPtr) storage_ = makePtr<VectorController<Real>>();
}

template<typename Real>
AffineTransformObjective<Real>::AffineTransformObjective(const Ptr<Objective<Real>>            &obj,
                                                         const Ptr<const LinearOperator<Real>> &A,
                                                         const Ptr<const Vector<Real>>         &b,
                                                         const Ptr<VectorController<Real>>     &storage)
  : obj_(obj), acon_(makePtr<LinearConstraint<Real>>(A,b)), storage_(storage) {
  primal_ = b->clone();
  Av_     = b->clone();
  dual_   = b->dual().clone();
  if (storage == nullPtr) storage_ = makePtr<VectorController<Real>>();
}

template<typename Real>
void AffineTransformObjective<Real>::update( const Vector<Real> &x, UpdateType type, int iter) {
  storage_->objectiveUpdate(type);
  acon_->update(x,type,iter);
  obj_->update(*transform(x),type,iter);
}

template<typename Real>
void AffineTransformObjective<Real>::update( const Vector<Real> &x, bool flag, int iter ) {
  storage_->objectiveUpdate(true);
  acon_->update(x,flag,iter);
  obj_->update(*transform(x),flag,iter);
}

template<typename Real>
Real AffineTransformObjective<Real>::value( const Vector<Real> &x, Real &tol ) {
  return obj_->value(*transform(x),tol); 
}

template<typename Real>
void AffineTransformObjective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  obj_->gradient(*dual_,*transform(x),tol);
  acon_->applyAdjointJacobian(g,*dual_,x,tol);
}

template<typename Real>
void AffineTransformObjective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  acon_->applyJacobian(*Av_,v,x,tol);
  obj_->hessVec(*dual_,*Av_,*transform(x),tol);
  acon_->applyAdjointJacobian(hv,*dual_,x,tol);
}

template<typename Real> 
void AffineTransformObjective<Real>::setParameter(const std::vector<Real> &param) {
  Objective<Real>::setParameter(param);
  obj_->setParameter(param);
}

template<typename Real>
Ptr<const Vector<Real>> AffineTransformObjective<Real>::transform(const Vector<Real> &x) {
  bool isApplied = storage_->get(*primal_,Objective<Real>::getParameter());
  if (!isApplied) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    acon_->value(*primal_,x,tol);
    storage_->set(*primal_,Objective<Real>::getParameter());
  }
  return primal_;
}

} // namespace ROL

#endif // ROL_AFFINE_TRANSFORM_OBJECTIVE_DEF_H
