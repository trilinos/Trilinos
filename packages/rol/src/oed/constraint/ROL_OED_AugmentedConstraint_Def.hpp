// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_AUGMENTED_PROBABILITY_CONSTRAINT_DEF_HPP
#define ROL_OED_AUGMENTED_PROBABILITY_CONSTRAINT_DEF_HPP

#include "ROL_ScalarLinearConstraint.hpp"
#include "ROL_OED_DesignVector.hpp"
#include "ROL_OED_ProbabilityConstraint.hpp"
#include "ROL_OED_ProbabilityConstraintRangeVector.hpp"

namespace ROL {
namespace OED {

template<typename Real>
AugmentedConstraint<Real>::AugmentedConstraint(const Vector<Real> &p,
                                               bool useScale,
                                               Real scale)
  : con_(makePtr<ProbabilityConstraint<Real>>(p,useScale,scale)),
    rvec_(makePtr<ProbabilityConstraintRangeVector<Real>>(
     dynamic_cast<const BatchStdVector<Real>&>(p).getBatchManager(),1)) {}
    //rvec_(makePtr<DualScaledStdVector<Real>>(
    //      makePtr<std::vector<Real>>(1,0),
    //      makePtr<std::vector<Real>>(1,static_cast<Real>(1)))) {}

template<typename Real>
AugmentedConstraint<Real>::AugmentedConstraint(const Ptr<Vector<Real>> &cost,
                                               Real budget)
  : con_(makePtr<ScalarLinearConstraint<Real>>(cost,budget)),
    rvec_(makePtr<SingletonVector<Real>>()) {}

template<typename Real>
void AugmentedConstraint<Real>::value(Vector<Real> &c,
           const Vector<Real> &x,
           Real &tol) {
  const DesignVector<Real> &xs = static_cast<const DesignVector<Real>&>(x);
  con_->value(c,*xs.getVector(),tol);
}

template<typename Real>
void AugmentedConstraint<Real>::applyJacobian(Vector<Real> &jv,
                   const Vector<Real> &v,
                   const Vector<Real> &x,
                   Real &tol) {
  const DesignVector<Real> &vs = static_cast<const DesignVector<Real>&>(v);
  const DesignVector<Real> &xs = static_cast<const DesignVector<Real>&>(x);
  con_->applyJacobian(jv,*vs.getVector(),*xs.getVector(),tol);
}

template<typename Real>
void AugmentedConstraint<Real>::applyAdjointJacobian(Vector<Real> &ajv,
                          const Vector<Real> &v,
                          const Vector<Real> &x,
                          Real &tol) {
  DesignVector<Real> &js = static_cast<DesignVector<Real>&>(ajv);
  const DesignVector<Real> &xs = static_cast<const DesignVector<Real>&>(x);
  con_->applyAdjointJacobian(*js.getVector(),v,*xs.getVector(),tol);
}

template<typename Real>
void AugmentedConstraint<Real>::applyAdjointHessian(Vector<Real> &ahuv,
                         const Vector<Real> &u,
                         const Vector<Real> &v,
                         const Vector<Real> &x,
                         Real &tol) {
  ahuv.zero();
}

} // End OED Namespace
} // End ROLNamespace

#endif
