// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_BOUND_CONSTRAINT_DEF_H
#define ROL_BOUND_CONSTRAINT_DEF_H

namespace ROL {

template<typename Real>
BoundConstraint<Real>::BoundConstraint(void)
  : Lactivated_(true), Uactivated_(true) {}

template<typename Real>
BoundConstraint<Real>::BoundConstraint(const Vector<Real> &x)
  : Lactivated_(false), Uactivated_(false) {
  try {
    lower_ = x.clone(); lower_->setScalar(ROL_NINF<Real>());
    upper_ = x.clone(); upper_->setScalar(ROL_INF<Real>());
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
    Ptr<Vector<Real>> Pv = v.clone();
    Pv->set(v);
    project(*Pv);
    Pv->axpy(static_cast<Real>(-1),v);
    Real diff = Pv->norm();
    return (diff <= ROL_EPSILON<Real>());
  }
  return true;
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
