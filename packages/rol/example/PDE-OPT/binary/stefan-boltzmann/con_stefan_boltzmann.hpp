// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  obj.hpp
    \brief Provides the interface for local (cell-based) objective function computations.
*/

#ifndef BINARY_CON_STEFAN_BOLTZMANN_HPP
#define BINARY_CON_STEFAN_BOLTZMANN_HPP

#include "ROL_StdConstraint.hpp"
#include "ROL_Bounds.hpp"

template <class Real>
class BudgetConstraint : public ROL::StdConstraint<Real> {
private:
  Real budget_;

public:
  BudgetConstraint(ROL::ParameterList &pl) {
    budget_ = pl.sublist("Problem").get("Control Budget",8.0);
  }

  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) override {
    c[0] = -budget_;
    for (const auto xv : x) c[0] += xv;
  }

  void applyJacobian( std::vector<Real> &jv, const std::vector<Real> &v, 
                      const std::vector<Real> &x, Real &tol ) override {
    jv[0] = static_cast<Real>(0);
    for (const auto vv : v) jv[0] += vv;
  }
  
  void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v, 
                             const std::vector<Real> &x, Real &tol ) override {
    ajv.assign(ajv.size(),v[0]);
  }

  void applyAdjointHessian( std::vector<Real> &ahuv, const std::vector<Real> &u,
                            const std::vector<Real> &v, const std::vector<Real> &x,
                            Real &tol ) {
    ahuv.assign(ahuv.size(),static_cast<Real>(0));
  }

  ROL::Ptr<ROL::Vector<Real>> createMultiplier(void) const {
    return ROL::makePtr<ROL::StdVector<Real>>(1,static_cast<Real>(0));
  }

  ROL::Ptr<ROL::BoundConstraint<Real>> createBounds(void) const {
    ROL::Ptr<ROL::Vector<Real>> l = createMultiplier(); l->setScalar(-budget_);
    ROL::Ptr<ROL::Vector<Real>> u = createMultiplier(); u->zero();
    return ROL::makePtr<ROL::Bounds<Real>>(l,u);
  }
}; // BudgetConstraint

#endif
