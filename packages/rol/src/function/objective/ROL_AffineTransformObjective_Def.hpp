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
  primal_ = acon_->createRangeSpaceVector();
  Av_     = acon_->createRangeSpaceVector();
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
