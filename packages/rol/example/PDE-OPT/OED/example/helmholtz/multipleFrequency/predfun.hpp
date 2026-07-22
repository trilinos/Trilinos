// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PDEOPT_PREDFUN_H
#define PDEOPT_PREDFUN_H

#include "model.hpp"

template <class Real>
class PredFun : public ROL::Objective<Real> {
private:
  const ROL::Ptr<ROL::Constraint<Real>> model_;
  const ROL::Ptr<ROL::Vector<Real>> vec_;
  const ROL::Ptr<ROL::Vector<Real>> vdual_;

public:
  PredFun(const ROL::Ptr<ROL::Constraint<Real>> &model,
          const ROL::Ptr<ROL::Vector<Real>> &vec)
    : model_(model), vec_(vec), vdual_(vec_->dual().clone()) {}

  void setParameter(const std::vector<Real> &param) {
    ROL::Objective<Real>::setParameter(param);
    model_->setParameter(param);
  }

  void update(const ROL::Vector<Real> &x, bool flag = true, int iter = -1) {
    model_->update(x,flag,iter);
  }

  Real value(const ROL::Vector<Real> &x, Real &tol ) {
    model_->value(*vdual_,x,tol);
    return vdual_->apply(*vec_);
  }

  void gradient(ROL::Vector<Real> &g,
                const ROL::Vector<Real> &x,
                Real &tol) {
    model_->applyAdjointJacobian(g,*vec_,x,tol);
  }

}; // class PredFun

#endif
