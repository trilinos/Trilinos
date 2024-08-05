// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CONSTRAINT_DYNAMICSTATE_H
#define ROL_CONSTRAINT_DYNAMICSTATE_H

#include "ROL_DynamicConstraint.hpp"
#include "ROL_Constraint.hpp"

namespace ROL {

template <class Real>
class Constraint_DynamicState : public Constraint<Real> {
private:
  const Ptr<DynamicConstraint<Real>> con_;
  const Ptr<const Vector<Real>>       uo_;
  const Ptr<const Vector<Real>>        z_;
  const Ptr<const TimeStamp<Real>>    ts_;

  Ptr<Vector<Real>> ijv_;
  bool isInit_;

public:
  Constraint_DynamicState(const Ptr<DynamicConstraint<Real>> &con,
                          const Ptr<const Vector<Real>>      &uo,
                          const Ptr<const Vector<Real>>      &z,
                          const Ptr<const TimeStamp<Real>>   &ts)
    : con_(con), uo_(uo), z_(z), ts_(ts), isInit_(false) {}

  void value(Vector<Real> &c,const Vector<Real> &u,Real &tol) {
    con_->value(c,*uo_,u,*z_,*ts_);
  }

  void applyJacobian(Vector<Real> &jv,const Vector<Real> &v,const Vector<Real> &u,Real &tol) {
    con_->applyJacobian_un(jv,v,*uo_,u,*z_,*ts_);
  }

  using Constraint<Real>::applyAdjointJacobian;
  void applyAdjointJacobian(Vector<Real> &ajv,const Vector<Real> &v,const Vector<Real> &u,Real &tol) {
    con_->applyAdjointJacobian_un(ajv,v,*uo_,u,*z_,*ts_);
  }

  void applyAdjointHessian(Vector<Real> &ahwv,const Vector<Real> &w,const Vector<Real> &v,const Vector<Real> &u,Real &tol) {
    con_->applyAdjointHessian_un_un(ahwv,w,v,*uo_,u,*z_,*ts_);
  }

  void update( const Vector<Real> &u, bool flag = true, int iter = -1 ) {
    //con_->update_un(u,*ts_);
    con_->update(*uo_,u,*z_,*ts_);
  }

  void applyPreconditioner(Vector<Real> &pv,const Vector<Real> &v,const Vector<Real> &u,const Vector<Real> &g,Real &tol) {
    if (!isInit_) {
      ijv_ = u.clone();
      isInit_ = true;
    }
    con_->applyInverseJacobian_un(*ijv_,v,*uo_,u,*z_,*ts_);
    con_->applyInverseAdjointJacobian_un(pv,ijv_->dual(),*uo_,u,*z_,*ts_);
  }

  // Definitions for parametrized (stochastic) equality constraints
  //void setParameter(const std::vector<Real> &param) {
  //  con_->setParameter(param);
  //  Constraint<Real>::setParameter(param);
  //}

}; // class Constraint_State

} // namespace ROL

#endif
