// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_RISKNEUTRALCONSTRAINT_HPP
#define ROL_RISKNEUTRALCONSTRAINT_HPP

#include "ROL_Ptr.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_SampleGenerator.hpp"

namespace ROL {

template<class Real>
class RiskNeutralConstraint : public Constraint<Real> {
private:
  const Ptr<Constraint<Real>>      con_;
  const Ptr<SampleGenerator<Real>> xsampler_;
  const Ptr<BatchManager<Real>>    cbman_;

  Ptr<Vector<Real>> conVec_;
  Ptr<Vector<Real>> optVec_;

  bool initialized_;

  void init(const Vector<Real> &c, const Vector<Real> &x) {
    if (!initialized_) {
      conVec_ = c.clone();
      optVec_ = x.dual().clone();
      initialized_ = true;
    }
  }

public:
  RiskNeutralConstraint( const Ptr<Constraint<Real>>      &con,
                         const Ptr<SampleGenerator<Real>> &xsampler,
                         const Ptr<BatchManager<Real>>    &cbman)
    : con_(con), xsampler_(xsampler), cbman_(cbman), initialized_(false) {}

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    con_->update(x,flag,iter);
  }

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) {
    con_->update(x,type,iter);
  }

  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol ) {
    init(c,x);
    conVec_->zero();
    for ( int i = 0; i < xsampler_->numMySamples(); ++i ) {
      con_->setParameter(xsampler_->getMyPoint(i));
      con_->value(c,x,tol);
      conVec_->axpy(xsampler_->getMyWeight(i),c);
    }
    c.zero();
    cbman_->sumAll(*conVec_,c);
  }

  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    init(jv,x);
    conVec_->zero();
    for ( int i = 0; i < xsampler_->numMySamples(); ++i ) {
      con_->setParameter(xsampler_->getMyPoint(i));
      con_->applyJacobian(jv,v,x,tol);
      conVec_->axpy(xsampler_->getMyWeight(i),jv);
    }
    jv.zero();
    cbman_->sumAll(*conVec_,jv);
  }

  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    init(v.dual(),x);
    optVec_->zero();
    for ( int i = 0; i < xsampler_->numMySamples(); ++i ) {
      con_->setParameter(xsampler_->getMyPoint(i));
      con_->applyAdjointJacobian(ajv,v,x,tol);
      optVec_->axpy(xsampler_->getMyWeight(i),ajv);
    }
    ajv.zero();
    xsampler_->sumAll(*optVec_,ajv);
  }

  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    init(u.dual(),x);
    optVec_->zero();
    for ( int i = 0; i < xsampler_->numMySamples(); ++i ) {
      con_->setParameter(xsampler_->getMyPoint(i));
      con_->applyAdjointHessian(ahuv,u,v,x,tol);
      optVec_->axpy(xsampler_->getMyWeight(i),ahuv);
    }
    ahuv.zero();
    xsampler_->sumAll(*optVec_,ahuv);
  }

};

}

#endif
