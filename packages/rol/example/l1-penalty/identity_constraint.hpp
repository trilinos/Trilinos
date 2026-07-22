// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef IDENTITY_CONSTRAINT_HPP
#define IDENTITY_CONSTRAINT_HPP

#include "ROL_Constraint.hpp"

namespace ROL {

template<typename Real>
class IdentityConstraint : public ROL::Constraint<Real> {
public:

  IdentityConstraint() = default;
  
  void value(       Vector<Real>& c,
              const Vector<Real>& x,
                    Real&         tol ) override {
    c.set(x);
  }

  void applyJacobian(       Vector<Real>& jv,
                      const Vector<Real>& v,
                      const Vector<Real>& x,
                            Real&         tol ) override {
    jv.set(v);
  }

  void applyAdjointJacobian(       Vector<Real>& ajv,
                             const Vector<Real>& v,
                             const Vector<Real>& x,
                                   Real&         tol ) override {
    ajv.set(v);
  }

  void applyAdjointHessian(       Vector<Real>& ahuv,
                            const Vector<Real>& u,
                            const Vector<Real>& v,
                            const Vector<Real>& x,
                                  Real&         tol ) override {
    ahuv.zero();
  }
};

} // namespace ROL

#endif //IDENTITY_CONSTRAINT_HPP

