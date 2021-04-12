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

#ifndef ROL_COMPOSITE_OBJECTIVE_DEF_H
#define ROL_COMPOSITE_OBJECTIVE_DEF_H

namespace ROL {

template<typename Real>
CompositeObjective<Real>::CompositeObjective(const std::vector<Ptr<Objective<Real>>> &obj_vec,
                                             const Ptr<StdObjective<Real>> &std_obj)
  : obj_vec_(obj_vec), std_obj_(std_obj), isInitialized_(false),
    isValueComputed_(false), isGradientComputed_(false) {
  obj_value_ = makePtr<std::vector<Real>>(obj_vec_.size(),0);
  obj_value_vec_ = makePtr<StdVector<Real>>(obj_value_);
  obj_grad_ = makePtr<std::vector<Real>>(obj_vec_.size(),0);
  obj_grad_vec_ = makePtr<StdVector<Real>>(obj_grad_);
  obj_gv_ = makePtr<std::vector<Real>>(obj_vec_.size(),0);
  obj_gv_vec_ = makePtr<StdVector<Real>>(obj_gv_);
  obj_hess_ = makePtr<std::vector<Real>>(obj_vec_.size(),0);
  obj_hess_vec_ = makePtr<StdVector<Real>>(obj_hess_);
}

template<typename Real>
void CompositeObjective<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {
  int size = obj_vec_.size();
  for (int i = 0; i < size; ++i) {
    obj_vec_[i]->update(x,type,iter);
  }
  isValueComputed_ = false;
  isGradientComputed_ = (type==UpdateType::Trial || type==UpdateType::Revert ? isGradientComputed_ : false);
}

template<typename Real>
void CompositeObjective<Real>::update( const Vector<Real> &x, bool flag, int iter ) {
  int size = obj_vec_.size();
  for (int i = 0; i < size; ++i) {
    obj_vec_[i]->update(x,flag,iter);
  }
  isValueComputed_ = false;
  isGradientComputed_ = (flag ? false : isGradientComputed_);
}

template<typename Real>
Real CompositeObjective<Real>::value( const Vector<Real> &x, Real &tol ) {
  computeValue(x,tol);
  return std_obj_->value(*obj_value_vec_,tol);
}

template<typename Real>
void CompositeObjective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  g.zero();
  computeGradient(x,tol);
  int size = obj_vec_.size();
  for (int i = 0; i < size; ++i) {
    g.axpy((*obj_grad_)[i],*(vec_grad_[i]));
  }
}

template<typename Real>
void CompositeObjective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  hv.zero();
  computeHessVec(v,x,tol);
  int size = obj_vec_.size();
  for (int i = 0; i < size; ++i) {
    hv.axpy((*obj_grad_)[i],*(vec_hess_[i]));
    hv.axpy((*obj_hess_)[i],*(vec_grad_[i]));
  }
}

template<typename Real>
void CompositeObjective<Real>::setParameter(const std::vector<Real> &param) {
  Objective<Real>::setParameter(param);
  const int size = obj_vec_.size();
  for (int i = 0; i < size; ++i) {
    obj_vec_[i]->setParameter(param);
  }
  std_obj_->setParameter(param);
  isValueComputed_ = false;    // Recompute value every time
  isGradientComputed_ = false; // Recompute gradient every time
}

template<typename Real>
void CompositeObjective<Real>::initialize(const Vector<Real> &x) {
  if (!isInitialized_){
    int size = obj_vec_.size();
    vec_grad_.clear(); vec_grad_.resize(size,nullPtr);
    vec_hess_.clear(); vec_hess_.resize(size,nullPtr);
    for (int i = 0; i < size; ++i) {
      vec_grad_[i] = x.dual().clone();
      vec_hess_[i] = x.dual().clone();
    }
    isInitialized_ = true;
  }
}

template<typename Real>
void CompositeObjective<Real>::computeValue(const Vector<Real> &x, Real &tol) {
  initialize(x);
  if (!isValueComputed_) {
    int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      (*obj_value_)[i] = obj_vec_[i]->value(x,tol);
    }
    isValueComputed_ = true;
  }
}

template<typename Real>
void CompositeObjective<Real>::computeGradient(const Vector<Real> &x, Real &tol) {
  computeValue(x,tol);
  if (!isGradientComputed_) {
    std_obj_->gradient(*(obj_grad_vec_),*(obj_value_vec_),tol);
    int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      obj_vec_[i]->gradient(*(vec_grad_[i]),x,tol);
    }
    isGradientComputed_ = true;
  }
}

template<typename Real>
void CompositeObjective<Real>::computeHessVec(const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
  computeGradient(x,tol);
  int size = obj_vec_.size();
  for (int i = 0; i < size; ++i) {
    //(*obj_gv_)[i] = vec_grad_[i]->dot(v.dual());
    (*obj_gv_)[i] = vec_grad_[i]->apply(v);
    obj_vec_[i]->hessVec(*(vec_hess_[i]),v,x,tol);
  }
  std_obj_->hessVec(*(obj_hess_vec_),*(obj_gv_vec_),*(obj_value_vec_),tol);
}

} // namespace ROL

#endif
