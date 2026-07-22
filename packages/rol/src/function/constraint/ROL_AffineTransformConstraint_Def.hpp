// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_AFFINE_TRANSFORM_CONSTRAINT_DEF_H
#define ROL_AFFINE_TRANSFORM_CONSTRAINT_DEF_H

namespace ROL {

template<typename Real>
AffineTransformConstraint<Real>::AffineTransformConstraint(const Ptr<Constraint<Real>>       &con,
                                                           const Ptr<Constraint<Real>>       &acon,
                                                           const Vector<Real>                &range,
                                                           const Ptr<VectorController<Real>> &storage)
  : con_(con), acon_(acon), storage_(storage) {
  primal_ = range.clone();
  Av_     = range.clone();
  dual_   = range.dual().clone();
  if (storage == nullPtr) storage_ = makePtr<VectorController<Real>>();
}

template<typename Real>
AffineTransformConstraint<Real>::AffineTransformConstraint(const Ptr<Constraint<Real>>       &con,
                                                           const Ptr<LinearConstraint<Real>> &acon,
                                                           const Ptr<VectorController<Real>> &storage)
  : con_(con), acon_(acon), storage_(storage) {
  primal_ = acon->createRangeSpaceVector();
  Av_     = acon->createRangeSpaceVector();
  dual_   = primal_->dual().clone();
  if (storage == nullPtr) storage_ = makePtr<VectorController<Real>>();
}

template<typename Real>
AffineTransformConstraint<Real>::AffineTransformConstraint(const Ptr<Constraint<Real>>           &con,
                                                           const Ptr<const LinearOperator<Real>> &A,
                                                           const Ptr<const Vector<Real>>         &b,
                                                           const Ptr<VectorController<Real>>     &storage)
  : con_(con), acon_(makePtr<LinearConstraint<Real>>(A,b)), storage_(storage) {
  primal_ = b->clone();
  Av_     = b->clone();
  dual_   = b->dual().clone();
  if (storage == nullPtr) storage_ = makePtr<VectorController<Real>>();
}

template<typename Real>
void AffineTransformConstraint<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {
  storage_->constraintUpdate(type);
  acon_->update(x,type,iter);
  con_->update(*transform(x),type,iter);
}

template<typename Real>
void AffineTransformConstraint<Real>::update( const Vector<Real> &x, bool flag, int iter ) {
  storage_->constraintUpdate(true);
  acon_->update(x,flag,iter);
  con_->update(*transform(x),flag,iter);
}

template<typename Real>
void AffineTransformConstraint<Real>::value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {
  con_->value(c,*transform(x),tol); 
}

template<typename Real>
void AffineTransformConstraint<Real>::applyJacobian( Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  acon_->applyJacobian(*Av_,v,x,tol);
  con_->applyJacobian(jv,*Av_,*transform(x),tol);
}

template<typename Real>
void AffineTransformConstraint<Real>::applyAdjointJacobian( Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  con_->applyAdjointJacobian(*dual_,v,*transform(x),tol);
  acon_->applyAdjointJacobian(ajv,*dual_,x,tol);
}

template<typename Real>
void AffineTransformConstraint<Real>::applyAdjointHessian( Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  acon_->applyJacobian(*Av_,v,x,tol);
  con_->applyAdjointHessian(*dual_,u,*Av_,*transform(x),tol);
  acon_->applyAdjointJacobian(ahuv,*dual_,x,tol);
}

template<typename Real>
Ptr<const Vector<Real>> AffineTransformConstraint<Real>::transform(const Vector<Real> &x) {
  bool isApplied = storage_->get(*primal_,Constraint<Real>::getParameter());
  if (!isApplied) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    acon_->value(*primal_,x,tol);
    storage_->set(*primal_,Constraint<Real>::getParameter());
  }
  return primal_;
}

} // namespace ROL

#endif // ROL_AFFINE_TRANSFORM_OBJECTIVE_H
