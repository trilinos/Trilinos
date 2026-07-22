// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_ALMOST_SURE_CONSTRAINT_H
#define ROL_ALMOST_SURE_CONSTRAINT_H

#include "ROL_RiskVector.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_SampleGenerator.hpp"
#include "ROL_SimulatedVector.hpp"

namespace ROL {

template <class Real>
class AlmostSureConstraint : public Constraint<Real> {
private:
  const Ptr<SampleGenerator<Real>> sampler_;
  const Ptr<Constraint<Real>>      con_;

  Ptr<Vector<Real>>                scratch1_;
  Ptr<Vector<Real>>                scratch2_;
  bool                             isInitialized_;

public:
  virtual ~AlmostSureConstraint() {}

  AlmostSureConstraint(const Ptr<SampleGenerator<Real>> &sampler,
                       const Ptr<Constraint<Real>> &con)
    : sampler_(sampler), con_(con), isInitialized_(false) {}

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    con_->update(x,flag,iter);
  }
  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) {
    con_->update(x,type,iter);
  }

  void value(Vector<Real> &c,
             const Vector<Real> &x,
             Real &tol) {
    c.zero();
    SimulatedVector<Real> &pc = dynamic_cast<SimulatedVector<Real>&>(c);

    std::vector<Real> param;
    for (int i = 0; i < sampler_->numMySamples(); ++i) {
      param  = sampler_->getMyPoint(i);
      con_->setParameter(param);
      con_->update(x);
      con_->value(*(pc.get(i)), x, tol);
    }
  }
 
  void applyJacobian(Vector<Real> &jv,
                     const Vector<Real> &v,
                     const Vector<Real> &x,
                     Real &tol) {
    jv.zero();
    SimulatedVector<Real> &pjv = dynamic_cast<SimulatedVector<Real>&>(jv);

    std::vector<Real> param;
    for (int i = 0; i < sampler_->numMySamples(); ++i) {
      param  = sampler_->getMyPoint(i);
      con_->setParameter(param);
      con_->update(x);
      con_->applyJacobian(*(pjv.get(i)), v, x, tol);
    }
  }

  void applyAdjointJacobian(Vector<Real> &ajv,
                            const Vector<Real> &v,
                            const Vector<Real> &x,
                            Real &tol) {
    ajv.zero();
    if (!isInitialized_) {
      scratch1_ = ajv.clone();
      scratch2_ = ajv.clone();
      isInitialized_ = true;
    }
    const SimulatedVector<Real> &pv = dynamic_cast<const SimulatedVector<Real>&>(v);

    std::vector<Real> param;
    scratch1_->zero(); scratch2_->zero();
    for (int i = 0; i < sampler_->numMySamples(); ++i) {
      param  = sampler_->getMyPoint(i);
      con_->setParameter(param);
      con_->update(x);
      con_->applyAdjointJacobian(*scratch1_, *(pv.get(i)), x, tol);
      scratch2_->plus(*scratch1_);
    }
    sampler_->sumAll(*scratch2_, ajv);
  }

  void applyAdjointHessian(Vector<Real> &ahuv,
                           const Vector<Real> &u,
                           const Vector<Real> &v,
                           const Vector<Real> &x,
                           Real &tol) {
    ahuv.zero();
    if (!isInitialized_) {
      scratch1_ = ahuv.clone();
      scratch2_ = ahuv.clone();
      isInitialized_ = true;
    }
    const SimulatedVector<Real> &pu = dynamic_cast<const SimulatedVector<Real>&>(u);

    std::vector<Real> param;
    scratch1_->zero(); scratch2_->zero();
    for (int i = 0; i < sampler_->numMySamples(); ++i) {
      param  = sampler_->getMyPoint(i);
      con_->setParameter(param);
      con_->update(x);
      con_->applyAdjointHessian(*scratch1_, *(pu.get(i)), v, x, tol);
      scratch2_->plus(*scratch1_);
    }
    sampler_->sumAll(*scratch2_, ahuv);
  }

  void applyPreconditioner(Vector<Real> &Pv,
                           const Vector<Real> &v,
                           const Vector<Real> &x,
                           const Vector<Real> &g,
                           Real &tol) {
    Pv.zero();
    SimulatedVector<Real> &ppv = dynamic_cast<SimulatedVector<Real>&>(Pv);
    const SimulatedVector<Real> &pv = dynamic_cast<const SimulatedVector<Real>&>(v);

    std::vector<Real> param;
    scratch1_->zero(); scratch2_->zero();
    for (int i = 0; i < sampler_->numMySamples(); ++i) {
      param  = sampler_->getMyPoint(i);
      con_->setParameter(param);
      con_->update(x);
      con_->applyPreconditioner(*(ppv.get(i)), *(pv.get(i)), x, g, tol);
    }
  }

}; // class AlmostSureConstraint

} // namespace ROL

#endif
