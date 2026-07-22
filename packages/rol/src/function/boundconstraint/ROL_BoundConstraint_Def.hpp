// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BOUND_CONSTRAINT_DEF_H
#define ROL_BOUND_CONSTRAINT_DEF_H

namespace ROL {

template<typename Real>
Real BoundConstraint<Real>::computeInf(const Vector<Real> &x) const {
  int dim = x.dimension();
  Real denom = (dim > 0 ? static_cast<Real>(dim) : 1e15);
  return std::sqrt(ROL_INF<Real>() / denom);
}

template<typename Real>
BoundConstraint<Real>::BoundConstraint(void)
  : Lactivated_(true), Uactivated_(true) {}

template<typename Real>
BoundConstraint<Real>::BoundConstraint(const Vector<Real> &x)
  : Lactivated_(false), Uactivated_(false) {
  try {
    lower_ = x.clone(); lower_->setScalar(-computeInf(x));
    upper_ = x.clone(); upper_->setScalar( computeInf(x));
  }
  catch(std::exception &e) {
    // Do nothing.  If someone calls getLowerBound or getUpperBound,
    // an exception will be thrown.
  }
}

template<typename Real>
void BoundConstraint<Real>::project( Vector<Real> &x ) {
  if (isActivated()) {
    throw Exception::NotImplemented(">>> ROL::BoundConstraint::project: Not Implemented!");
  }
}

template<typename Real>
void BoundConstraint<Real>::projectInterior( Vector<Real> &x ) {
  if (isActivated()) {
    throw Exception::NotImplemented(">>> ROL::BoundConstraint::projectInterior: Not Implemented!");
  }
}

template<typename Real>
void BoundConstraint<Real>::pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps ) {
  if (isUpperActivated()) {
    throw Exception::NotImplemented(">>> ROL::BoundConstraint::pruneUpperActive: Not Implemented!");
  }
}

template<typename Real>
void BoundConstraint<Real>::pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps, Real geps ) {
  if (isUpperActivated()) {
    throw Exception::NotImplemented(">>> ROL::BoundConstraint::pruneUpperActive: Not Implemented!");
  }
}

template<typename Real>
void BoundConstraint<Real>::pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps ) {
  if (isLowerActivated()) {
    throw Exception::NotImplemented(">>> ROL::BoundConstraint::pruneLowerActive: Not Implemented!");
  }
}

template<typename Real>
void BoundConstraint<Real>::pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps, Real geps ) {
  if (isLowerActivated()) {
    throw Exception::NotImplemented(">>> ROL::BoundConstraint::pruneLowerActive: Not Implemented!");
  }
}

template<typename Real>
const Ptr<const Vector<Real>> BoundConstraint<Real>::getLowerBound( void ) const {
  if (lower_ != nullPtr) {
    return lower_;
  }
  throw Exception::NotImplemented(">>> ROL::BoundConstraint::getLowerBound: Lower bound not provided!");
}

template<typename Real>
const Ptr<const Vector<Real>> BoundConstraint<Real>::getUpperBound( void ) const {
  if (upper_ != nullPtr) {
    return upper_;
  }
  throw Exception::NotImplemented(">>> ROL::BoundConstraint::getUpperBound: Upper bound not provided!");
}

template<typename Real>
bool BoundConstraint<Real>::isFeasible( const Vector<Real> &v ) { 
  if (isActivated()) {
    const Real tol(static_cast<Real>(1e-2)*std::sqrt(ROL_EPSILON<Real>()));
    Ptr<Vector<Real>> Pv = v.clone();
    Pv->set(v);
    project(*Pv);
    Pv->axpy(static_cast<Real>(-1),v);
    Real diff = Pv->norm();
    return (diff <= tol);
  }
  return true;
}

template<typename Real>
void BoundConstraint<Real>::applyInverseScalingFunction(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const {
  throw Exception::NotImplemented(">>> BoundConstraint::applyInverseScalingFunction : This function has not been implemeted!");
}

template<typename Real>
void BoundConstraint<Real>::applyScalingFunctionJacobian(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const {
  throw Exception::NotImplemented(">>> BoundConstraint::applyScalingFunctionJacobian : This function has not been implemeted!");
}

template<typename Real>
void BoundConstraint<Real>::activateLower(void) {
  Lactivated_ = true;
}

template<typename Real>
void BoundConstraint<Real>::activateUpper(void) {
  Uactivated_ = true;
}

template<typename Real>
void BoundConstraint<Real>::activate(void) {
  activateLower();
  activateUpper();
}

template<typename Real>
void BoundConstraint<Real>::deactivateLower(void) {
  Lactivated_ = false;
}

template<typename Real>
void BoundConstraint<Real>::deactivateUpper(void) {
  Uactivated_ = false;
}

template<typename Real>
void BoundConstraint<Real>::deactivate(void) {
  deactivateLower();
  deactivateUpper();
}

template<typename Real>
bool BoundConstraint<Real>::isLowerActivated(void) const {
  return Lactivated_;
}

template<typename Real>
bool BoundConstraint<Real>::isUpperActivated(void) const {
  return Uactivated_;
}

template<typename Real>
bool BoundConstraint<Real>::isActivated(void) const {
  return (isLowerActivated() || isUpperActivated());
}

template<typename Real>
void BoundConstraint<Real>::pruneActive( Vector<Real> &v, const Vector<Real> &x, Real eps ) {
  if (isActivated()) {
    pruneUpperActive(v,x,eps);
    pruneLowerActive(v,x,eps);
  }
}

template<typename Real>
void BoundConstraint<Real>::pruneActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps, Real geps ) {
  if (isActivated()) {
    pruneUpperActive(v,g,x,xeps,geps);
    pruneLowerActive(v,g,x,xeps,geps);
  }
}

template<typename Real>
void BoundConstraint<Real>::pruneLowerInactive( Vector<Real> &v, const Vector<Real> &x, Real eps ) {
  if (isLowerActivated()) {
    const Real one(1);
    Ptr<Vector<Real>> tmp = v.clone(); 
    tmp->set(v);
    pruneLowerActive(*tmp,x,eps);
    v.axpy(-one,*tmp);
  }
}

template<typename Real>
void BoundConstraint<Real>::pruneUpperInactive( Vector<Real> &v, const Vector<Real> &x, Real eps ) { 
  if (isUpperActivated()) {
    const Real one(1);
    Ptr<Vector<Real>> tmp = v.clone(); 
    tmp->set(v);
    pruneUpperActive(*tmp,x,eps);
    v.axpy(-one,*tmp);
  }
}

template<typename Real>
void BoundConstraint<Real>::pruneLowerInactive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps, Real geps ) { 
  if (isLowerActivated()) {
    const Real one(1);
    Ptr<Vector<Real>> tmp = v.clone(); 
    tmp->set(v);
    pruneLowerActive(*tmp,g,x,xeps,geps);
    v.axpy(-one,*tmp);
  }
}

template<typename Real>
void BoundConstraint<Real>::pruneUpperInactive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps, Real geps ) { 
  if (isUpperActivated()) {
    const Real one(1);
    Ptr<Vector<Real>> tmp = v.clone(); 
    tmp->set(v);
    pruneUpperActive(*tmp,g,x,xeps,geps);
    v.axpy(-one,*tmp);
  }
}

template<typename Real>
void BoundConstraint<Real>::pruneInactive( Vector<Real> &v, const Vector<Real> &x, Real eps ) { 
  if (isActivated()) {
    const Real one(1);
    Ptr<Vector<Real>> tmp = v.clone(); 
    tmp->set(v);
    pruneActive(*tmp,x,eps);
    v.axpy(-one,*tmp);
  }
}

template<typename Real>
void BoundConstraint<Real>::pruneInactive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps, Real geps ) { 
  if (isActivated()) {
    const Real one(1);
    Ptr<Vector<Real>> tmp = v.clone(); 
    tmp->set(v);
    pruneActive(*tmp,g,x,xeps,geps);
    v.axpy(-one,*tmp);
  }
}

template<typename Real>
void BoundConstraint<Real>::computeProjectedGradient( Vector<Real> &g, const Vector<Real> &x ) {
  if (isActivated()) {
    Ptr<Vector<Real>> tmp = g.clone();
    tmp->set(g);
    pruneActive(g,*tmp,x);
  }
}

template<typename Real>
void BoundConstraint<Real>::computeProjectedStep( Vector<Real> &v, const Vector<Real> &x ) { 
  if (isActivated()) {
    const Real one(1);
    v.plus(x);
    project(v);
    v.axpy(-one,x);
  }
}

} // namespace ROL

#endif
