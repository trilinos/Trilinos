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
