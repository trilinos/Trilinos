// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BINARY_CONSTRAINT_DEF_H
#define ROL_BINARY_CONSTRAINT_DEF_H

namespace ROL {

template<typename Real>
BinaryConstraint<Real>::BinaryConstraint( const ROL::Ptr<const Vector<Real>> &lo,
                  const ROL::Ptr<const Vector<Real>> &up, Real gamma ) :
    lo_(lo), up_(up), d_(lo_->clone()), gamma_(gamma) {} 

template<typename Real>
BinaryConstraint<Real>::BinaryConstraint( const BoundConstraint<Real> &bnd, Real gamma ) :
    BinaryConstraint( bnd.getLowerBound(), bnd.getUpperBound(), gamma ) {}
 
template<typename Real>
BinaryConstraint<Real>::BinaryConstraint( const ROL::Ptr<const BoundConstraint<Real>> &bnd, Real gamma ) :
    BinaryConstraint( bnd->getLowerBound(), bnd->getUpperBound(), gamma ) {}

template<typename Real>
void BinaryConstraint<Real>::value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
  const Real one(1);
  c.set( x );
  c.axpy( -one, *lo_ ); // c = x-l
  d_->set( *up_ );
  d_->axpy( -one, x );  // d = u-x
  c.applyBinary(BoundsCheck(0), *d_ );
  c.scale( gamma_ );
}

template<typename Real>
void BinaryConstraint<Real>::applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
  const Real one(1);
  jv.set( x );
  jv.axpy( -one, *lo_ );
  d_->set( *up_ );
  d_->axpy( -one, x );
  jv.applyBinary( BoundsCheck(1), *d_ );
  jv.applyBinary( Elementwise::Multiply<Real>(), v );
  jv.scale( gamma_ );
}

template<typename Real>
void BinaryConstraint<Real>::applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
  applyJacobian(ajv,v,x,tol); 
}

template<typename Real>
void BinaryConstraint<Real>::applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
  const Real one(1);
  ahuv.set( x );
  ahuv.axpy( -one, *lo_ );
  d_->set( *up_ );
  d_->axpy( -one, x );
  ahuv.applyBinary( BoundsCheck(2), *d_ );
  ahuv.applyBinary( Elementwise::Multiply<Real>(), v );
  ahuv.applyBinary( Elementwise::Multiply<Real>(), u );
  ahuv.scale( gamma_ ); 
}

template<typename Real>
void BinaryConstraint<Real>::setPenalty( Real gamma ) {
  gamma_ = gamma;
}

} // namespace ROL

#endif // ROL_BINARY_CONSTRAINT_DEF_H
