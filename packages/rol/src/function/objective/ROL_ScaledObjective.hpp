// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SCALED_OBJECTIVE_HPP
#define ROL_SCALED_OBJECTIVE_HPP

#include "ROL_Objective.hpp"

namespace ROL {

template <typename Real>
class ScaledObjective : public Objective<Real> {
private:
  const Ptr<Objective<Real>> obj_;
  const Real scale_;

public:
  ScaledObjective(const Ptr<Objective<Real>> &obj, Real scale)
    : obj_(obj), scale_(scale) {}

  void update(const Vector<Real> &x, UpdateType type, int iter = -1) override;
  void setParameter(const std::vector<Real> &param) override;
  Real value( const Vector<Real> &x, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;
  void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;
  void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;

}; // class ScaledObjective

} // End ROL Namespace

#include "ROL_ScaledObjective_Def.hpp"

#endif
