// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_BINARYCONSTRAINT_DEF_HPP
#define ROL_OED_BINARYCONSTRAINT_DEF_HPP

namespace ROL::OED {

template<typename Real>
void BinaryConstraint<Real>::value(
       std::vector<Real> &c,const std::vector<Real> &x,Real &tol) {
  startTimer("value");
  const Real one(1);
  const unsigned dim=x.size();
  for(unsigned i=0u; i<dim; ++i) c[i] = x[i]*(x[i]-one);
  stopTimer("value");
}

template<typename Real>
void BinaryConstraint<Real>::applyJacobian(
       std::vector<Real> &jv,const std::vector<Real> &v,
       const std::vector<Real> &x,Real &tol) {
  startTimer("applyJacobian");
  const Real one(1), two(2);
  const unsigned dim=x.size();
  for(unsigned i=0u; i<dim; ++i) jv[i] = (two*x[i]-one)*v[i];
  stopTimer("applyJacobian");
}

template<typename Real>
void BinaryConstraint<Real>::applyAdjointJacobian(
       std::vector<Real> &ajv,const std::vector<Real> &v,
       const std::vector<Real> &x,Real &tol) {
  startTimer("applyAdjointJacobian");
  const Real one(1), two(2);
  const unsigned dim=x.size();
  for(unsigned i=0u; i<dim; ++i) ajv[i] = (two*x[i]-one)*v[i];
  stopTimer("applyAdjointJacobian");
}

template<typename Real>
void BinaryConstraint<Real>::applyAdjointHessian(
       std::vector<Real> &ahwv,const std::vector<Real> &w,
       const std::vector<Real> &v,const std::vector<Real> &x,Real &tol) {
  startTimer("applyAdjointHessian");
  const Real two(2);
  const unsigned dim=x.size();
  for(unsigned i=0u; i<dim; ++i) ahwv[i] = two*w[i]*v[i];
  stopTimer("applyAdjointHessian");
}

} // End ROL::OED Namespace

#endif
