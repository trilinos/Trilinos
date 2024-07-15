// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PDEOPT_MODEL_H
#define PDEOPT_MODEL_H

#include "ROL_StdObjective.hpp"
#include "ROL_SampledScalar.hpp"
#include "ROL_SampledVector.hpp"

template <class Real>
class Model : public ROL::StdObjective<Real> {
private:
  const ROL::Ptr<ROL::Objective<Real>> obs_;
  const std::vector<ROL::Ptr<ROL::Vector<Real>>> state_;

  std::vector<ROL::Ptr<ROL::SampledScalar<Real>>> value_storage_;
  ROL::Ptr<ROL::SampledVector<Real>> grad_storage_;
  ROL::Ptr<ROL::Vector<Real>> g_;

  Real computeValue(const int index) {
    Real tol = std::sqrt(ROL::ROL_EPSILON<Real>()), val(0);
    std::vector<Real> param = ROL::StdObjective<Real>::getParameter();
    bool isComputed = value_storage_[index]->get(val,param);
    if (!isComputed) {
      bool isComputedG = grad_storage_->get(*g_,param);
      if (!isComputedG) {
        obs_->setParameter(param);
        obs_->gradient(*g_,*state_[index],tol);
        grad_storage_->set(*g_,param);
      }
      val = state_[index]->dot(g_->dual());
      value_storage_[index]->set(val,param);
    }
    return val;
  }

public:
  Model(const ROL::Ptr<ROL::Objective<Real>> &obs,
        const std::vector<ROL::Ptr<ROL::Vector<Real>>> &state)
    : obs_(obs), state_(state) {
    int nfactors = state_.size();
    value_storage_.clear(); value_storage_.resize(nfactors);
    for (int i = 0; i < nfactors; ++i) {
      value_storage_[i] = ROL::makePtr<ROL::SampledScalar<Real>>();
    }
    grad_storage_ = ROL::makePtr<ROL::SampledVector<Real>>();
    g_ = state[0]->dual().clone();
  }

  void setParameter(const std::vector<Real> &param) {
    ROL::StdObjective<Real>::setParameter(param);
    obs_->setParameter(param);
  }

  Real value( const std::vector<Real> &x, Real &tol ) {
    int nfactors = state_.size();
    Real val(0);
    for (int i = 0; i < nfactors; ++i) {
      val += computeValue(i) * x[i];
    }
    return val;
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    int nfactors = state_.size();
    for (int i = 0; i < nfactors; ++i) {
      g[i] = computeValue(i);
    }
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real &tol ) {
    hv.zero();
  }

}; // class Model

#endif
