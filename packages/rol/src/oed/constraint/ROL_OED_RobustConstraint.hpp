// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_ROBUSTCONSTRAINT_HPP
#define ROL_OED_ROBUSTCONSTRAINT_HPP

#include "ROL_Ptr.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_ScaledStdVector.hpp"
#include "ROL_OED_ObjectiveArray.hpp"
#include "ROL_OED_DesignVector.hpp"

namespace ROL {
namespace OED {

template<typename Real>
class RobustConstraint : public Constraint<Real> {
private:
  const Ptr<ObjectiveArray<Real>> objVec_;
  const Ptr<Vector<Real>> dwa_;

public:
  RobustConstraint(const Ptr<ObjectiveArray<Real>> &objVec)
    : objVec_(objVec), dwa_(((objVec->buildDesignVector())->dual()).clone()) {}

  const Ptr<DesignVector<Real>> buildDomainVector() const {
    return makePtr<DesignVector<Real>>(objVec_->buildDesignVector());
  }

  const Ptr<StdVector<Real>> buildRangeVector() const {
    auto data = makePtr<std::vector<Real>>(objVec_->numObjectives(),static_cast<Real>(0));
    auto scal = makePtr<std::vector<Real>>(objVec_->getWeights());
    return makePtr<PrimalScaledStdVector<Real>>(data,scal);
  }

  void update(const Vector<Real> &x, UpdateType type, int iter = -1 ) override {
    const DesignVector<Real> &xs = static_cast<const DesignVector<Real>&>(x);
    for (unsigned i = 0u; i < objVec_->numObjectives(); ++i)
      objVec_->getObjective(i)->update(*xs.getVector(),type,iter);
  }

  void value(Vector<Real> &c,const Vector<Real> &x,Real &tol) override {
    StdVector<Real> &cs = static_cast<StdVector<Real>&>(c);
    const DesignVector<Real> &xs = static_cast<const DesignVector<Real>&>(x);
    for (unsigned i = 0u; i < objVec_->numObjectives(); ++i)
      (*cs.getVector())[i] = objVec_->getObjective(i)->value(*xs.getVector(),tol) - xs.getValue();
  }

  void applyJacobian(Vector<Real> &jv,const Vector<Real> &v,const Vector<Real> &x,Real &tol) override {
    StdVector<Real> &js = static_cast<StdVector<Real>&>(jv);
    const DesignVector<Real> &vs = static_cast<const DesignVector<Real>&>(v);
    const DesignVector<Real> &xs = static_cast<const DesignVector<Real>&>(x);
    for (unsigned i = 0u; i < objVec_->numObjectives(); ++i) {
      objVec_->getObjective(i)->gradient(*dwa_,*xs.getVector(),tol);
      (*js.getVector())[i] = dwa_->apply(*vs.getVector()) - vs.getValue();
    }
  }

  void applyAdjointJacobian(Vector<Real> &ajv,const Vector<Real> &v,const Vector<Real> &x,Real &tol) override {
    ajv.zero();
    DesignVector<Real> &js = static_cast<DesignVector<Real>&>(ajv);
    const StdVector<Real> &vs = static_cast<const StdVector<Real>&>(v);
    const DesignVector<Real> &xs = static_cast<const DesignVector<Real>&>(x);
    Real val(0);
    for (unsigned i = 0u; i < objVec_->numObjectives(); ++i) {
      objVec_->getObjective(i)->gradient(*dwa_,*xs.getVector(),tol);
      js.getVector()->axpy((*vs.getVector())[i],*dwa_);
      val += (*vs.getVector())[i];
    }
    js.setValue(-val);
  }

  void applyAdjointHessian(Vector<Real> &ahwv,const Vector<Real> &w, const Vector<Real> &v,
                           const Vector<Real> &x,Real &tol) override {
    ahwv.zero();
    DesignVector<Real> &hs = static_cast<DesignVector<Real>&>(ahwv);
    const StdVector<Real> &ws = static_cast<const StdVector<Real>&>(w);
    const DesignVector<Real> &vs = static_cast<const DesignVector<Real>&>(v);
    const DesignVector<Real> &xs = static_cast<const DesignVector<Real>&>(x);
    for (unsigned i = 0u; i < objVec_->numObjectives(); ++i) {
      objVec_->getObjective(i)->hessVec(*dwa_,*vs.getVector(),*xs.getVector(),tol);
      hs.getVector()->axpy((*ws.getVector())[i],*dwa_);
    }
  }

  void setParameter(const std::vector<Real> &param) {
    Constraint<Real>::setParameter(param);
    for (unsigned i = 0u; i < objVec_->numObjectives(); ++i)
      objVec_->getObjective(i)->setParameter(param);
  }

}; // class RobustConstraint

} // End OED Namespace
} // End ROL Namespace

//#include "ROL_OED_RobustConstraint_Def.hpp"

#endif
