// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SLACKLESSOBJECTIVE_DEF_HPP
#define ROL_SLACKLESSOBJECTIVE_DEF_HPP

namespace ROL {

template<typename Real> 
SlacklessObjective<Real>::SlacklessObjective( const Ptr<Objective<Real>> &obj ) : obj_(obj) {}

template<typename Real>
Ptr<Objective<Real>> SlacklessObjective<Real>::getObjective(void) const {
  return obj_;
}

template<typename Real> 
void SlacklessObjective<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {
  obj_->update( *getOpt(x), type, iter );
}

template<typename Real> 
void SlacklessObjective<Real>::update( const Vector<Real> &x, bool flag, int iter ) {
  obj_->update( *getOpt(x), flag, iter );
}

template<typename Real> 
Real SlacklessObjective<Real>::value( const Vector<Real> &x, Real &tol ) {
  return obj_->value( *getOpt(x), tol );
}

template<typename Real> 
Real SlacklessObjective<Real>::dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol ) {
  return obj_->dirDeriv(*getOpt(x),*getOpt(d),tol);
}

template<typename Real> 
void SlacklessObjective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  zeroSlack(g);
  obj_->gradient(*getOpt(g),*getOpt(x),tol);
}

template<typename Real> 
void SlacklessObjective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  zeroSlack(hv);
  obj_->hessVec(*getOpt(hv),*getOpt(v),*getOpt(x),tol);     
}

template<typename Real> 
void SlacklessObjective<Real>::invHessVec( Vector<Real> &ihv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  obj_->invHessVec( *getOpt(ihv), *getOpt(v), *getOpt(x), tol );
  PartitionedVector<Real> &Pvp = dynamic_cast<PartitionedVector<Real>&>(ihv);
  const PartitionedVector<Real> &vp = dynamic_cast<const PartitionedVector<Real>&>(v);
  const int nvec = static_cast<int>(vp.numVectors());
  for (int i = 1; i < nvec; ++i) {
    Pvp.get(i)->set(vp.get(i)->dual());
  }
}

template<typename Real> 
void SlacklessObjective<Real>::precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  obj_->precond( *getOpt(Pv), *getOpt(v), *getOpt(x), tol );
  PartitionedVector<Real> &Pvp = dynamic_cast<PartitionedVector<Real>&>(Pv);
  const PartitionedVector<Real> &vp = dynamic_cast<const PartitionedVector<Real>&>(v);
  const int nvec = static_cast<int>(vp.numVectors());
  for (int i = 1; i < nvec; ++i) {
    Pvp.get(i)->set(vp.get(i)->dual());
  }
}

template<typename Real> 
void SlacklessObjective<Real>::setParameter(const std::vector<Real> &param) {
  Objective<Real>::setParameter(param);
  obj_->setParameter(param);
}

template<typename Real> 
Ptr<Vector<Real>> SlacklessObjective<Real>::getOpt( Vector<Real> &xs ) const {
  return dynamic_cast<PartitionedVector<Real>&>(xs).get(0);
}

template<typename Real> 
Ptr<const Vector<Real>> SlacklessObjective<Real>::getOpt( const Vector<Real> &xs ) const {
  return dynamic_cast<const PartitionedVector<Real>&>(xs).get(0);
}

template<typename Real> 
void SlacklessObjective<Real>::zeroSlack( Vector<Real> &x ) const {
  PartitionedVector<Real> &xpv
    = dynamic_cast<PartitionedVector<Real>&>(x);
  const int nvec = static_cast<int>(xpv.numVectors());
  for (int i = 1; i < nvec; ++i) {
    xpv.get(i)->zero();
  }
} 

} // namespace ROL

#endif // ROL__SLACKLESSOBJECTIVE_HPP

