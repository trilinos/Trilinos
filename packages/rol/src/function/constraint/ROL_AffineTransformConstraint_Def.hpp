// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
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
