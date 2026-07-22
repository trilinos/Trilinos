// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SLACKLESSCONSTRAINT_DEF_HPP
#define ROL_SLACKLESSCONSTRAINT_DEF_HPP

namespace ROL {

template<typename Real> 
SlacklessConstraint<Real>::SlacklessConstraint( const Ptr<Constraint<Real>> &con ) : con_(con) {}
 
template<typename Real> 
void SlacklessConstraint<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {
  con_->update( *getOpt(x), type, iter );
}
 
template<typename Real> 
void SlacklessConstraint<Real>::update( const Vector<Real> &x, bool flag, int iter ) {
  con_->update( *getOpt(x), flag, iter );
}

template<typename Real> 
void SlacklessConstraint<Real>::value(Vector<Real> &c, const Vector<Real> &x, Real &tol ) {
  con_->value( c, *getOpt(x), tol );
}

template<typename Real> 
void SlacklessConstraint<Real>::applyJacobian( Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  con_->applyJacobian( jv, *getOpt(v), *getOpt(x), tol );
}

template<typename Real> 
void SlacklessConstraint<Real>::applyAdjointJacobian( Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &dualv, Real &tol ) {
  zeroSlack(ajv);
  con_->applyAdjointJacobian( *getOpt(ajv), v, *getOpt(x), dualv, tol );
}

template<typename Real> 
void SlacklessConstraint<Real>::applyAdjointHessian( Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  zeroSlack(ahuv);
  con_->applyAdjointHessian( *getOpt(ahuv), u, *getOpt(v), *getOpt(x), tol );     
}

template<typename Real> 
void SlacklessConstraint<Real>::setParameter(const std::vector<Real> &param) {
  Constraint<Real>::setParameter(param);
  con_->setParameter(param);
}

template<typename Real> 
Ptr<Vector<Real>> SlacklessConstraint<Real>::getOpt( Vector<Real> &xs ) const {
  return dynamic_cast<PartitionedVector<Real>&>(xs).get(0);
}

template<typename Real> 
Ptr<const Vector<Real>> SlacklessConstraint<Real>::getOpt( const Vector<Real> &xs ) const {
  return dynamic_cast<const PartitionedVector<Real>&>(xs).get(0);
}

template<typename Real> 
void SlacklessConstraint<Real>::zeroSlack( Vector<Real> &x ) const {
  PartitionedVector<Real> &xpv
    = dynamic_cast<PartitionedVector<Real>&>(x);
  const int nvec = static_cast<int>(xpv.numVectors());
  for (int i = 1; i < nvec; ++i) {
    xpv.get(i)->zero();
  }
} 

} // namespace ROL

#endif // ROL__SLACKLESSCONSTRAINT_HPP

