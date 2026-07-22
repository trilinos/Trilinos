// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_COMPOSITEOBJECTIVE_SIMOPT_H
#define ROL_COMPOSITEOBJECTIVE_SIMOPT_H

#include "ROL_StdObjective.hpp"
#include "ROL_Objective_SimOpt.hpp"

/** @ingroup func_group
    \class ROL::CompositeObjective_SimOpt
    \brief Provides the interface to evaluate simulation-based composite
    objective functions.
*/

namespace ROL {

template <class Real>
class CompositeObjective_SimOpt : public Objective_SimOpt<Real> {
private:
  const std::vector<ROL::Ptr<Objective_SimOpt<Real> > > obj_vec_;
  const ROL::Ptr<StdObjective<Real> > std_obj_;

  ROL::Ptr<std::vector<Real>> obj_value_;
  ROL::Ptr<std::vector<Real>> obj_grad_;
  ROL::Ptr<std::vector<Real>> obj_gv_;
  ROL::Ptr<std::vector<Real>> obj_hess_;
  ROL::Ptr<StdVector<Real>> obj_value_vec_;
  ROL::Ptr<StdVector<Real>> obj_grad_vec_;
  ROL::Ptr<StdVector<Real>> obj_gv_vec_;
  ROL::Ptr<StdVector<Real>> obj_hess_vec_;
  std::vector<ROL::Ptr<Vector<Real>>> vec_grad1_;
  std::vector<ROL::Ptr<Vector<Real>>> vec_grad2_;
  std::vector<ROL::Ptr<Vector<Real>>> vec_hess1_;
  std::vector<ROL::Ptr<Vector<Real>>> vec_hess2_;

  bool isInitialized_, isValueComputed_;
  bool isGradientComputed_, isGradient1Computed_, isGradient2Computed_;

  void initialize(const Vector<Real> &u, const Vector<Real> &z) {
    if (!isInitialized_){
      int size = obj_vec_.size();
      vec_grad1_.clear(); vec_grad1_.resize(size,ROL::nullPtr);
      vec_grad2_.clear(); vec_grad2_.resize(size,ROL::nullPtr);
      vec_hess1_.clear(); vec_hess1_.resize(size,ROL::nullPtr);
      vec_hess2_.clear(); vec_hess2_.resize(size,ROL::nullPtr);
      for (int i = 0; i < size; ++i) {
        vec_grad1_[i] = u.dual().clone();
        vec_grad2_[i] = z.dual().clone();
        vec_hess1_[i] = u.dual().clone();
        vec_hess2_[i] = z.dual().clone();
      }
      isInitialized_ = true;
    }
  }

  void computeValue(const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
    initialize(u,z);
    if (!isValueComputed_) {
      int size = obj_vec_.size();
      for (int i = 0; i < size; ++i) {
        (*obj_value_)[i] = obj_vec_[i]->value(u,z,tol);
      }
      isValueComputed_ = true;
    }
  }

  void computeGradient(const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
    computeValue(u,z,tol);
    if (!isGradientComputed_) {
      std_obj_->gradient(*(obj_grad_vec_),*(obj_value_vec_),tol);
      isGradientComputed_ = true;
    }
  }

  void computeGradient1(const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
    computeGradient(u,z,tol);
    if (!isGradient1Computed_) {
      int size = obj_vec_.size();
      for (int i = 0; i < size; ++i) {
        obj_vec_[i]->gradient_1(*(vec_grad1_[i]),u,z,tol);
      }
      isGradient1Computed_ = true;
    }
  }

  void computeGradient2(const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
    computeGradient(u,z,tol);
    if (!isGradient2Computed_) {
      int size = obj_vec_.size();
      for (int i = 0; i < size; ++i) {
        obj_vec_[i]->gradient_2(*(vec_grad2_[i]),u,z,tol);
      }
      isGradient2Computed_ = true;
    }
  }

  void computeHessVec11(const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
    computeGradient1(u,z,tol);
    int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      //(*obj_gv_)[i] = vec_grad1_[i]->dot(v.dual());
      (*obj_gv_)[i] = vec_grad1_[i]->apply(v);
      obj_vec_[i]->hessVec_11(*(vec_hess1_[i]),v,u,z,tol);
    }
    std_obj_->hessVec(*(obj_hess_vec_),*(obj_gv_vec_),*(obj_value_vec_),tol);
  }

  void computeHessVec12(const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
    computeGradient1(u,z,tol);
    computeGradient2(u,z,tol);
    int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      //(*obj_gv_)[i] = vec_grad2_[i]->dot(v.dual());
      (*obj_gv_)[i] = vec_grad2_[i]->apply(v);
      obj_vec_[i]->hessVec_12(*(vec_hess1_[i]),v,u,z,tol);
    }
    std_obj_->hessVec(*(obj_hess_vec_),*(obj_gv_vec_),*(obj_value_vec_),tol);
  }

  void computeHessVec21(const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
    computeGradient1(u,z,tol);
    computeGradient2(u,z,tol);
    int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      //(*obj_gv_)[i] = vec_grad1_[i]->dot(v.dual());
      (*obj_gv_)[i] = vec_grad1_[i]->apply(v);
      obj_vec_[i]->hessVec_21(*(vec_hess2_[i]),v,u,z,tol);
    }
    std_obj_->hessVec(*(obj_hess_vec_),*(obj_gv_vec_),*(obj_value_vec_),tol);
  }

  void computeHessVec22(const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
    computeGradient2(u,z,tol);
    int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      //(*obj_gv_)[i] = vec_grad2_[i]->dot(v.dual());
      (*obj_gv_)[i] = vec_grad2_[i]->apply(v);
      obj_vec_[i]->hessVec_22(*(vec_hess2_[i]),v,u,z,tol);
    }
    std_obj_->hessVec(*(obj_hess_vec_),*(obj_gv_vec_),*(obj_value_vec_),tol);
  }

public:
  CompositeObjective_SimOpt(const std::vector<ROL::Ptr<Objective_SimOpt<Real> > > &obj_vec,
                            const ROL::Ptr<StdObjective<Real> > &std_obj)
    : obj_vec_(obj_vec), std_obj_(std_obj),
      isInitialized_(false), isValueComputed_(false),
      isGradientComputed_(false), isGradient1Computed_(false), isGradient2Computed_(false) {
    obj_value_ = ROL::makePtr<std::vector<Real>>(obj_vec_.size(),0);
    obj_value_vec_ = ROL::makePtr<StdVector<Real>>(obj_value_);
    obj_grad_ = ROL::makePtr<std::vector<Real>>(obj_vec_.size(),0);
    obj_grad_vec_ = ROL::makePtr<StdVector<Real>>(obj_grad_);
    obj_gv_ = ROL::makePtr<std::vector<Real>>(obj_vec_.size(),0);
    obj_gv_vec_ = ROL::makePtr<StdVector<Real>>(obj_gv_);
    obj_hess_ = ROL::makePtr<std::vector<Real>>(obj_vec_.size(),0);
    obj_hess_vec_ = ROL::makePtr<StdVector<Real>>(obj_hess_);
  }

  void update( const Vector<Real> &u, const Vector<Real> &z, bool flag = true, int iter = -1 ) {
    int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      obj_vec_[i]->update(u,z,flag,iter);
    }
    isValueComputed_     = false;
    isGradientComputed_  = (flag ? false : isGradientComputed_);
    isGradient1Computed_ = (flag ? false : isGradient1Computed_);
    isGradient2Computed_ = (flag ? false : isGradient2Computed_);
  }

  void update( const Vector<Real> &u, const Vector<Real> &z, UpdateType type, int iter = -1 ) {
    int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      obj_vec_[i]->update(u,z,type,iter);
    }
    // Do something smarter here
    isValueComputed_     = false;
    isGradientComputed_  = false;
    isGradient1Computed_ = false;
    isGradient2Computed_ = false;
  }

  Real value( const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    computeValue(u,z,tol);
    return std_obj_->value(*obj_value_vec_,tol);
  }


  void gradient_1( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    g.zero();
    computeGradient1(u,z,tol);
    int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      g.axpy((*obj_grad_)[i],*(vec_grad1_[i]));
    }
  }

  void gradient_2( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    g.zero();
    computeGradient2(u,z,tol);
    int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      g.axpy((*obj_grad_)[i],*(vec_grad2_[i]));
    }
  }

  void hessVec_11( Vector<Real> &hv, const Vector<Real> &v, 
             const Vector<Real> &u,  const Vector<Real> &z, Real &tol ) {
    hv.zero();
    computeHessVec11(v,u,z,tol);
    int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      hv.axpy((*obj_grad_)[i],*(vec_hess1_[i]));
      hv.axpy((*obj_hess_)[i],*(vec_grad1_[i]));
    }
  }

  void hessVec_12( Vector<Real> &hv, const Vector<Real> &v, 
             const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    hv.zero();
    computeHessVec12(v,u,z,tol);
    int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      hv.axpy((*obj_grad_)[i],*(vec_hess1_[i]));
      hv.axpy((*obj_hess_)[i],*(vec_grad1_[i]));
    }
  }

  void hessVec_21( Vector<Real> &hv, const Vector<Real> &v, 
             const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    hv.zero();
    computeHessVec21(v,u,z,tol);
    int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      hv.axpy((*obj_grad_)[i],*(vec_hess2_[i]));
      hv.axpy((*obj_hess_)[i],*(vec_grad2_[i]));
    }
  }

  void hessVec_22( Vector<Real> &hv, const Vector<Real> &v, 
             const Vector<Real> &u,  const Vector<Real> &z, Real &tol ) {
    hv.zero();
    computeHessVec22(v,u,z,tol);
    int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      hv.axpy((*obj_grad_)[i],*(vec_hess2_[i]));
      hv.axpy((*obj_hess_)[i],*(vec_grad2_[i]));
    }
  }

// Definitions for parametrized (stochastic) objective functions
public:
  void setParameter(const std::vector<Real> &param) {
    Objective_SimOpt<Real>::setParameter(param);
    const int size = obj_vec_.size();
    for (int i = 0; i < size; ++i) {
      obj_vec_[i]->setParameter(param);
    }
    std_obj_->setParameter(param);
    isValueComputed_ = false;    // Recompute value every time
    isGradientComputed_ = false; // Recompute gradient every time
    isGradient1Computed_ = false; // Recompute gradient every time
    isGradient2Computed_ = false; // Recompute gradient every time
  }
};

} // namespace ROL

#endif
