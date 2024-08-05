// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STOCHASTIC_CONSTRAINT_H
#define ROL_STOCHASTIC_CONSTRAINT_H

#include "ROL_StochasticObjective.hpp"
#include "ROL_ConstraintFromObjective.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template <class Real>
class StochasticConstraint : public Constraint<Real> {
private:
  Ptr<StochasticObjective<Real>> robj_;
  Ptr<Constraint<Real>>          con_;
  Ptr<SampleGenerator<Real>>     sampler_;

public:
  StochasticConstraint(const Ptr<Objective<Real>> &obj,
               const Ptr<SampleGenerator<Real>>   &sampler,
               ParameterList                      &parlist,
               const int index = 0)
    : sampler_(sampler) {
    robj_ = makePtr<StochasticObjective<Real>>(obj,parlist,sampler,1,index);
    con_  = makePtr<ConstraintFromObjective<Real>>(robj_);
  }

  StochasticConstraint(const Ptr<Constraint<Real>> &con,
               const Ptr<SampleGenerator<Real>>    &sampler,
               ParameterList                       &parlist,
               const int index = 0)
    : sampler_(sampler) {
    try {
      Ptr<ConstraintFromObjective<Real>> cfo
        = dynamicPtrCast<ConstraintFromObjective<Real>>(con);
      robj_ = makePtr<StochasticObjective<Real>>(cfo->getObjective(),
                parlist,sampler,1,index);
      con_  = makePtr<ConstraintFromObjective<Real>>(robj_);
    }
    catch (std::exception &e) {
      throw Exception::NotImplemented(">>> ROL::StochasticConstraint: Input constraint must be a ConstraintFromObjective!");
    }
  }

  Real computeStatistic(const Vector<Real> &x) const {
    return robj_->computeStatistic(x);
  }

  void setIndex(int ind) {
    robj_->setIndex(ind);
  }

  void update(const Vector<Real> &x, UpdateType type, int iter = -1) {
    con_->update(x,type,iter);
  }

  void update(const Vector<Real> &x, bool flag = true, int iter = -1) {
    con_->update(x,flag,iter);
  }

  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
    con_->value(c,x,tol);
  }

  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    con_->applyJacobian(jv,v,x,tol);
  }

  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    con_->applyAdjointJacobian(ajv,v,x,tol);
  }

  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    con_->applyAdjointHessian(ahuv,u,v,x,tol);
  }

}; // class StochasticConstraint

} // namespace ROL

#endif
