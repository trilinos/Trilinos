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

#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"

template <class Real>
class Model : public ROL::Constraint<Real> {
private:
  const std::vector<ROL::Ptr<ROL::Objective<Real>>> objs_;
  ROL::Ptr<ROL::Vector<Real>> xdual_;

public:
  Model(const std::vector<ROL::Ptr<ROL::Objective<Real>>> &objs)
    : objs_(objs) {}

  void setParameter(const std::vector<Real> &param) {
    ROL::Constraint<Real>::setParameter(param);
    for (const auto obj : objs_) obj->setParameter(param);
  }

  void update(const ROL::Vector<Real> &x, bool flag = true, int iter = -1) {
    for (const auto obj : objs_) obj->update(x,flag,iter);
  }

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &x, Real &tol ) {
    ROL::Ptr<std::vector<Real>> cdata = static_cast<ROL::StdVector<Real>&>(c).getVector();
    int nobs = objs_.size();
    for (int i = 0; i < nobs; ++i) {
      (*cdata)[i] = objs_[i]->value(x,tol);
    }
  }

  void applyJacobian(ROL::Vector<Real> &jv,
                     const ROL::Vector<Real> &v,
                     const ROL::Vector<Real> &x,
                     Real &tol) {
    if (xdual_==ROL::nullPtr) xdual_ = x.dual().clone();
    ROL::Ptr<std::vector<Real>> jdata = static_cast<ROL::StdVector<Real>&>(jv).getVector();
    int nobs = objs_.size();
    for (int i = 0; i < nobs; ++i) {
      xdual_->zero();
      objs_[i]->gradient(*xdual_,x,tol);
      (*jdata)[i] = xdual_->dot(v.dual());
    }
  }

  void applyAdjointJacobian(ROL::Vector<Real> &ajv,
                            const ROL::Vector<Real> &v,
                            const ROL::Vector<Real> &x,
                            Real &tol) {
    if (xdual_==ROL::nullPtr) xdual_ = x.dual().clone();
    ROL::Ptr<const std::vector<Real>> vdata = static_cast<const ROL::StdVector<Real>&>(v).getVector();
    ajv.zero();
    int nobs = objs_.size();
    for (int i = 0; i < nobs; ++i) {
      xdual_->zero();
      objs_[i]->gradient(*xdual_,x,tol);
      ajv.axpy((*vdata)[i],*xdual_);
    }
  }

}; // class Model

#endif
