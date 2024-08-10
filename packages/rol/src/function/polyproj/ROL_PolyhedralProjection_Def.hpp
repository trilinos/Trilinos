// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_POLYHEDRALPROJECTION_DEF_H
#define ROL_POLYHEDRALPROJECTION_DEF_H

namespace ROL {

template<typename Real>
PolyhedralProjection<Real>::PolyhedralProjection(const Ptr<BoundConstraint<Real>> &bnd)
  : bnd_(bnd), con_(nullPtr) {}

template<typename Real>
PolyhedralProjection<Real>::PolyhedralProjection(const Vector<Real>               &xprim,
                                                 const Vector<Real>               &xdual,
                                                 const Ptr<BoundConstraint<Real>> &bnd,
                                                 const Ptr<Constraint<Real>>      &con,
                                                 const Vector<Real>               &mul,
                                                 const Vector<Real>               &res)
  : bnd_(bnd), con_(con) {
  xprim_ = xprim.clone();
  xdual_ = xdual.clone();
  mul_   = mul.clone();
  res_   = res.clone();
}

template<typename Real>
void PolyhedralProjection<Real>::project(Vector<Real> &x, std::ostream &stream) {
  if (con_ == nullPtr) {
    bnd_->project(x);
  }
  else {
    throw Exception::NotImplemented(">>> ROL::PolyhedralProjection::project : No projection implemented!");
  }
}

template<typename Real>
const Ptr<Constraint<Real>> PolyhedralProjection<Real>::getLinearConstraint(void) const {
  return con_;
}

template<typename Real>
const Ptr<BoundConstraint<Real>> PolyhedralProjection<Real>::getBoundConstraint(void) const {
  return bnd_;
}

template<typename Real>
const Ptr<Vector<Real>> PolyhedralProjection<Real>::getMultiplier(void) const {
  return mul_;
}

template<typename Real>
const Ptr<Vector<Real>> PolyhedralProjection<Real>::getResidual(void) const {
  return res_;
}

} // namespace ROL

#endif
