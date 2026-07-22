// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_RISKLESSOBJECTIVE_DEF_HPP
#define ROL_RISKLESSOBJECTIVE_DEF_HPP

namespace ROL {

template<typename Real>
RiskLessObjective<Real>::RiskLessObjective(const Ptr<Objective<Real>> &obj) : obj_(obj) {}

template<typename Real>
void RiskLessObjective<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {
  Ptr<const Vector<Real>> x0
    = dynamic_cast<const RiskVector<Real>&>(x).getVector();
  obj_->update(*x0,type,iter);
}

template<typename Real>
void RiskLessObjective<Real>::update( const Vector<Real> &x, bool flag, int iter ) {
  Ptr<const Vector<Real>> x0
    = dynamic_cast<const RiskVector<Real>&>(x).getVector();
  obj_->update(*x0,flag,iter);
}

template<typename Real>
Real RiskLessObjective<Real>::value( const Vector<Real> &x, Real &tol ) {
  Ptr<const Vector<Real>> x0
    = dynamic_cast<const RiskVector<Real>&>(x).getVector();
  return obj_->value(*x0,tol);
}

template<typename Real>
void RiskLessObjective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  Ptr<Vector<Real>> g0
    = dynamic_cast<RiskVector<Real>&>(g).getVector();
  Ptr<const Vector<Real>> x0
    = dynamic_cast<const RiskVector<Real>&>(x).getVector();
  obj_->gradient(*g0,*x0,tol);
}

template<typename Real>
void RiskLessObjective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v,
              const Vector<Real> &x, Real &tol ) {
  Ptr<Vector<Real>> hv0
    = dynamic_cast<RiskVector<Real>&>(hv).getVector();
  Ptr<const Vector<Real>> v0
    = dynamic_cast<const RiskVector<Real>&>(v).getVector();
  Ptr<const Vector<Real>> x0
    = dynamic_cast<const RiskVector<Real>&>(x).getVector();
  obj_->hessVec(*hv0,*v0,*x0,tol);
}

template<typename Real>
void RiskLessObjective<Real>::precond( Vector<Real> &Pv, const Vector<Real> &v,
              const Vector<Real> &x, Real &tol ) {
  Ptr<Vector<Real>> Pv0
    = dynamic_cast<RiskVector<Real>&>(Pv).getVector();
  Ptr<const Vector<Real>> v0
    = dynamic_cast<const RiskVector<Real>&>(v).getVector();
  Ptr<const Vector<Real>> x0
    = dynamic_cast<const RiskVector<Real>&>(x).getVector();
  obj_->precond(*Pv0,*v0,*x0,tol);
}

template<typename Real>
void RiskLessObjective<Real>::setParameter(const std::vector<Real> &param) {
  Objective<Real>::setParameter(param);
  obj_->setParameter(param);
}

}

#endif
