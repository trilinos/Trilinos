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

#ifndef ROL_COMPOSITE_OBJECTIVE_H
#define ROL_COMPOSITE_OBJECTIVE_H

#include "ROL_StdObjective.hpp"

/** @ingroup func_group
    \class ROL::CompositeObjective
    \brief Provides the interface to evaluate composite objective functions.
*/


namespace ROL {

template <class Real>
class CompositeObjective : public Objective<Real> {
private:
  const std::vector<ROL::Ptr<Objective<Real> > > obj_vec_;
  const ROL::Ptr<StdObjective<Real> > std_obj_;

  ROL::Ptr<std::vector<Real> > obj_value_, obj_grad_, obj_gv_, obj_hess_;
  ROL::Ptr<StdVector<Real> > obj_value_vec_, obj_grad_vec_, obj_gv_vec_, obj_hess_vec_;
  std::vector<ROL::Ptr<Vector<Real> > > vec_grad_, vec_hess_;

  bool isInitialized_, isValueComputed_, isGradientComputed_;

  void initialize(const Vector<Real> &x) {
    if (!isInitialized_){
      int size = obj_vec_.size();
      vec_grad_.clear(); vec_grad_.resize(size,ROL::nullPtr);
      vec_hess_.clear(); vec_hess_.resize(size,ROL::nullPtr);
      for (int i = 0; i < size; ++i) {
        vec_grad_[i] = x.dual().clone();
        vec_hess_[i] = x.dual().clone();
      }
      isInitialized_ = true;
    }
  }

  void computeValue(const Vector<Real> &x, Real &tol) {
    initialize(x);
    if (!isValueComputed_) {
      int size = obj_vec_.size();
      for (int i = 0; i < size; ++i) {
        (*obj_value_)[i] = obj_vec_[i]->value(x,tol);
      }
      isValueComputed_ = true;
    }
  }

  void computeGradient(const Vector<Real> &x, Real &tol) {
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

  void computeHessVec(const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    computeGradient(x,tol);
    int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      (*obj_gv_)[i] = vec_grad_[i]->dot(v.dual());
      obj_vec_[i]->hessVec(*(vec_hess_[i]),v,x,tol);
    }
    std_obj_->hessVec(*(obj_hess_vec_),*(obj_gv_vec_),*(obj_value_vec_),tol);
  }

public:
  CompositeObjective(const std::vector<ROL::Ptr<Objective<Real> > > &obj_vec,
                     const ROL::Ptr<StdObjective<Real> > &std_obj)
    : obj_vec_(obj_vec), std_obj_(std_obj), isInitialized_(false),
      isValueComputed_(false), isGradientComputed_(false) {
    obj_value_ = ROL::makePtr<std::vector<Real>>(obj_vec_.size(),0);
    obj_value_vec_ = ROL::makePtr<StdVector<Real>>(obj_value_);
    obj_grad_ = ROL::makePtr<std::vector<Real>>(obj_vec_.size(),0);
    obj_grad_vec_ = ROL::makePtr<StdVector<Real>>(obj_grad_);
    obj_gv_ = ROL::makePtr<std::vector<Real>>(obj_vec_.size(),0);
    obj_gv_vec_ = ROL::makePtr<StdVector<Real>>(obj_gv_);
    obj_hess_ = ROL::makePtr<std::vector<Real>>(obj_vec_.size(),0);
    obj_hess_vec_ = ROL::makePtr<StdVector<Real>>(obj_hess_);
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      obj_vec_[i]->update(x,flag,iter);
    }
    isValueComputed_ = false;
    isGradientComputed_ = (flag ? false : isGradientComputed_);
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    computeValue(x,tol);
    return std_obj_->value(*obj_value_vec_,tol);
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    g.zero();
    computeGradient(x,tol);
    int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      g.axpy((*obj_grad_)[i],*(vec_grad_[i]));
    }
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    hv.zero();
    computeHessVec(v,x,tol);
    int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      hv.axpy((*obj_grad_)[i],*(vec_hess_[i]));
      hv.axpy((*obj_hess_)[i],*(vec_grad_[i]));
    }
  }

// Definitions for parametrized (stochastic) objective functions
public:
  void setParameter(const std::vector<Real> &param) {
    Objective<Real>::setParameter(param);
    const int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      obj_vec_[i]->setParameter(param);
    }
    std_obj_->setParameter(param);
    isValueComputed_ = false;    // Recompute value every time
    isGradientComputed_ = false; // Recompute gradient every time
  }
};

} // namespace ROL

#endif
