// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SIMULATED_CONSTRAINT_H
#define ROL_SIMULATED_CONSTRAINT_H

#include "ROL_SimulatedVector.hpp"
#include "ROL_RiskVector.hpp"
#include "ROL_Constraint_SimOpt.hpp"

namespace ROL {

template <class Real>
class SimulatedConstraint : public Constraint<Real> {
private:
  const ROL::Ptr<SampleGenerator<Real> > sampler_;
  const ROL::Ptr<Constraint_SimOpt<Real> > pcon_;
  const bool useWeights_;

public:

  virtual ~SimulatedConstraint() {}

  SimulatedConstraint(const ROL::Ptr<SampleGenerator<Real> > & sampler,
                              const ROL::Ptr<Constraint_SimOpt<Real> > & pcon,
                              const bool useWeights = true)
    : sampler_(sampler), pcon_(pcon), useWeights_(useWeights) {}

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    pcon_->update(x,flag,iter);
  }
  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) {
    pcon_->update(x,type,iter);
  }

  void value(Vector<Real> &c,
             const Vector<Real> &x,
             Real &tol) {
    c.zero();
    SimulatedVector<Real> &pc = dynamic_cast<SimulatedVector<Real>&>(c);
    const Vector_SimOpt<Real> &uz = dynamic_cast<const Vector_SimOpt<Real>&>(x);
    ROL::Ptr<const Vector<Real> > uptr = uz.get_1();
    ROL::Ptr<const Vector<Real> > zptr = uz.get_2();
    try {
      const RiskVector<Real> &rz = dynamic_cast<const RiskVector<Real>&>(*zptr);
      zptr = rz.getVector();
    }
    catch (const std::bad_cast &e) {}
    const SimulatedVector<Real> &pu = dynamic_cast<const SimulatedVector<Real>&>(*uptr);

    std::vector<Real> param;
    Real weight(0), one(1);
    for (typename std::vector<SimulatedVector<Real> >::size_type i=0; i<pu.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pcon_->setParameter(param);
      pcon_->update(*(pu.get(i)), *zptr);
      pcon_->value(*(pc.get(i)), *(pu.get(i)), *zptr, tol);
      weight = (useWeights_) ? weight : one;
      pc.get(i)->scale(weight);
    }

  }

 
  virtual void applyJacobian(Vector<Real> &jv,
                             const Vector<Real> &v,
                             const Vector<Real> &x,
                             Real &tol) {
    jv.zero();
    // cast jv
    SimulatedVector<Real> &pjv = dynamic_cast<SimulatedVector<Real>&>(jv);
    // split x
    const Vector_SimOpt<Real> &xuz = dynamic_cast<const Vector_SimOpt<Real>&>(x);
    ROL::Ptr<const Vector<Real> > xuptr = xuz.get_1();
    ROL::Ptr<const Vector<Real> > xzptr = xuz.get_2();
    try {
      const RiskVector<Real> &rxz = dynamic_cast<const RiskVector<Real>&>(*xzptr);
      xzptr = rxz.getVector();
    }
    catch (const std::bad_cast &e) {}
    const SimulatedVector<Real> &pxu = dynamic_cast<const SimulatedVector<Real>&>(*xuptr);
    // split v
    const Vector_SimOpt<Real> &vuz = dynamic_cast<const Vector_SimOpt<Real>&>(v);
    ROL::Ptr<const Vector<Real> > vuptr = vuz.get_1();
    ROL::Ptr<const Vector<Real> > vzptr = vuz.get_2();
    try {
      const RiskVector<Real> &rvz = dynamic_cast<const RiskVector<Real>&>(*vzptr);
      vzptr = rvz.getVector();
    }
    catch (const std::bad_cast &e) {}
    const SimulatedVector<Real> &pvu = dynamic_cast<const SimulatedVector<Real>&>(*vuptr);

    std::vector<Real> param;
    Real weight(0), one(1);
    for (typename std::vector<SimulatedVector<Real> >::size_type i=0; i<pvu.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pcon_->setParameter(param);
      Vector_SimOpt<Real> xi(ROL::constPtrCast<Vector<Real> >(pxu.get(i)), ROL::constPtrCast<Vector<Real> >(xzptr));
      Vector_SimOpt<Real> vi(ROL::constPtrCast<Vector<Real> >(pvu.get(i)), ROL::constPtrCast<Vector<Real> >(vzptr));
      pcon_->update(xi);
      pcon_->applyJacobian(*(pjv.get(i)), vi, xi, tol);
      weight = (useWeights_) ? weight : one;
      pjv.get(i)->scale(weight);
    }
  }


  virtual void applyAdjointJacobian(Vector<Real> &ajv,
                                    const Vector<Real> &v,
                                    const Vector<Real> &x,
                                    Real &tol) {
    ajv.zero();
    // split ajv
    Vector_SimOpt<Real> &ajvuz = dynamic_cast<Vector_SimOpt<Real>&>(ajv);
    ROL::Ptr<Vector<Real> > ajvuptr = ajvuz.get_1();
    ROL::Ptr<Vector<Real> > ajvzptr = ajvuz.get_2();
    try {
      RiskVector<Real> &rajvz = dynamic_cast<RiskVector<Real>&>(*ajvzptr);
      ajvzptr = rajvz.getVector();
    }
    catch (const std::bad_cast &e) {}
    SimulatedVector<Real> &pajvu = dynamic_cast<SimulatedVector<Real>&>(*ajvuptr);
    // split x
    const Vector_SimOpt<Real> &xuz = dynamic_cast<const Vector_SimOpt<Real>&>(x);
    ROL::Ptr<const Vector<Real> > xuptr = xuz.get_1();
    ROL::Ptr<const Vector<Real> > xzptr = xuz.get_2();
    try {
      const RiskVector<Real> &rxz = dynamic_cast<const RiskVector<Real>&>(*xzptr);
      xzptr = rxz.getVector();
    }
    catch (const std::bad_cast &e) {}
    const SimulatedVector<Real> &pxu = dynamic_cast<const SimulatedVector<Real>&>(*xuptr);
    // cast v
    const SimulatedVector<Real> &pv = dynamic_cast<const SimulatedVector<Real>&>(v);

    std::vector<Real> param;
    Real weight(0), one(1);
    ROL::Ptr<Vector<Real> > tmp1 = ajvzptr->clone();
    ROL::Ptr<Vector<Real> > tmp2 = ajvzptr->clone();
    for (typename std::vector<SimulatedVector<Real> >::size_type i=0; i<pv.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pcon_->setParameter(param);
      Vector_SimOpt<Real> xi(ROL::constPtrCast<Vector<Real> >(pxu.get(i)), ROL::constPtrCast<Vector<Real> >(xzptr));
      Vector_SimOpt<Real> ajvi(pajvu.get(i), tmp1);
      pcon_->update(xi);
      pcon_->applyAdjointJacobian(ajvi, *(pv.get(i)), xi, tol);
      weight = (useWeights_) ? weight : one;
      ajvi.scale(weight);
      tmp2->plus(*tmp1);
    }
    sampler_->sumAll(*tmp2, *ajvzptr);

  }


  virtual void applyAdjointHessian(Vector<Real> &ahuv,
                                   const Vector<Real> &u,
                                   const Vector<Real> &v,
                                   const Vector<Real> &x,
                                   Real &tol) {
    ahuv.zero();
    // split ahuv
    Vector_SimOpt<Real> &ahuvuz = dynamic_cast<Vector_SimOpt<Real>&>(ahuv);
    ROL::Ptr<Vector<Real> > ahuvuptr = ahuvuz.get_1();
    ROL::Ptr<Vector<Real> > ahuvzptr = ahuvuz.get_2();
    try {
      RiskVector<Real> &rahuvz = dynamic_cast<RiskVector<Real>&>(*ahuvzptr);
      ahuvzptr = rahuvz.getVector();
    }
    catch (const std::bad_cast &e) {}
    SimulatedVector<Real> &pahuvu = dynamic_cast<SimulatedVector<Real>&>(*ahuvuptr);
    // cast u
    const SimulatedVector<Real> &pu = dynamic_cast<const SimulatedVector<Real>&>(u);
    // split v
    const Vector_SimOpt<Real> &vuz = dynamic_cast<const Vector_SimOpt<Real>&>(v);
    ROL::Ptr<const Vector<Real> > vuptr = vuz.get_1();
    ROL::Ptr<const Vector<Real> > vzptr = vuz.get_2();
    try {
      const RiskVector<Real> &rvz = dynamic_cast<const RiskVector<Real>&>(*vzptr);
      vzptr = rvz.getVector();
    }
    catch (const std::bad_cast &e) {}
    const SimulatedVector<Real> &pvu = dynamic_cast<const SimulatedVector<Real>&>(*vuptr);
    // split x
    const Vector_SimOpt<Real> &xuz = dynamic_cast<const Vector_SimOpt<Real>&>(x);
    ROL::Ptr<const Vector<Real> > xuptr = xuz.get_1();
    ROL::Ptr<const Vector<Real> > xzptr = xuz.get_2();
    try {
      const RiskVector<Real> &rxz = dynamic_cast<const RiskVector<Real>&>(*xzptr);
      xzptr = rxz.getVector();
    }
    catch (const std::bad_cast &e) {}
    const SimulatedVector<Real> &pxu = dynamic_cast<const SimulatedVector<Real>&>(*xuptr);

    std::vector<Real> param;
    Real weight(0), one(1);
    ROL::Ptr<Vector<Real> > tmp1 = ahuvzptr->clone();
    ROL::Ptr<Vector<Real> > tmp2 = ahuvzptr->clone();
    for (typename std::vector<SimulatedVector<Real> >::size_type i=0; i<pxu.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pcon_->setParameter(param);
      Vector_SimOpt<Real> xi(ROL::constPtrCast<Vector<Real> >(pxu.get(i)), ROL::constPtrCast<Vector<Real> >(xzptr));
      Vector_SimOpt<Real> vi(ROL::constPtrCast<Vector<Real> >(pvu.get(i)), ROL::constPtrCast<Vector<Real> >(vzptr));
      Vector_SimOpt<Real> ahuvi(pahuvu.get(i), tmp1);
      pcon_->update(xi);
      pcon_->applyAdjointHessian(ahuvi, *(pu.get(i)), vi, xi, tol);
      weight = (useWeights_) ? weight : one;
      ahuvi.scale(weight);
      tmp2->plus(*tmp1);
    }
    sampler_->sumAll(*tmp2, *ahuvzptr);

  }

  virtual void applyPreconditioner(Vector<Real> &Pv,
                                   const Vector<Real> &v,
                                   const Vector<Real> &x,
                                   const Vector<Real> &g,
                                   Real &tol) {
    Pv.zero();
    // cast Pv
    SimulatedVector<Real> &ppv = dynamic_cast<SimulatedVector<Real>&>(Pv);
    // split x
    const Vector_SimOpt<Real> &xuz = dynamic_cast<const Vector_SimOpt<Real>&>(x);
    ROL::Ptr<const Vector<Real> > xuptr = xuz.get_1();
    ROL::Ptr<const Vector<Real> > xzptr = xuz.get_2();
    try {
      const RiskVector<Real> &rxz = dynamic_cast<const RiskVector<Real>&>(*xzptr);
      xzptr = rxz.getVector();
    }
    catch (const std::bad_cast &e) {}
    const SimulatedVector<Real> &pxu = dynamic_cast<const SimulatedVector<Real>&>(*xuptr);
    // split g
    const Vector_SimOpt<Real> &guz = dynamic_cast<const Vector_SimOpt<Real>&>(g);
    ROL::Ptr<const Vector<Real> > guptr = guz.get_1();
    ROL::Ptr<const Vector<Real> > gzptr = guz.get_2();
    try {
      const RiskVector<Real> &rgz = dynamic_cast<const RiskVector<Real>&>(*gzptr);
      gzptr = rgz.getVector();
    }
    catch (const std::bad_cast &e) {}
    const SimulatedVector<Real> &pgu = dynamic_cast<const SimulatedVector<Real>&>(*guptr);
    // cast v
    const SimulatedVector<Real> &pv = dynamic_cast<const SimulatedVector<Real>&>(v);

    std::vector<Real> param;
    Real weight(0), one(1);
    for (typename std::vector<SimulatedVector<Real> >::size_type i=0; i<pv.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pcon_->setParameter(param);
      Vector_SimOpt<Real> xi(ROL::constPtrCast<Vector<Real> >(pxu.get(i)), ROL::constPtrCast<Vector<Real> >(xzptr));
      Vector_SimOpt<Real> gi(ROL::constPtrCast<Vector<Real> >(pgu.get(i)), ROL::constPtrCast<Vector<Real> >(gzptr));
      pcon_->update(xi);
      pcon_->applyPreconditioner(*(ppv.get(i)), *(pv.get(i)), xi, gi, tol);
      weight = (useWeights_) ? weight : one;
      ppv.get(i)->scale(one/(weight*weight));
    }

  }


}; // class SimulatedConstraint

} // namespace ROL

#endif
