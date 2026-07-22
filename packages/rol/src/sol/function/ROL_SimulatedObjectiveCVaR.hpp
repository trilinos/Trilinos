// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SIMULATED_OBJECTIVE_CVAR_H
#define ROL_SIMULATED_OBJECTIVE_CVAR_H

#include "ROL_SimulatedVector.hpp"
#include "ROL_PlusFunction.hpp"
#include "ROL_RiskVector.hpp"
#include "ROL_Objective_SimOpt.hpp"

namespace ROL {

template <class Real>
class SimulatedObjectiveCVaR : public Objective<Real> {
private:
  const ROL::Ptr<SampleGenerator<Real> > sampler_;
  const ROL::Ptr<Objective_SimOpt<Real> > pobj_;
  const ROL::Ptr<PlusFunction<Real> > pfunc_;
  const Real alpha_;

public:

  virtual ~SimulatedObjectiveCVaR() {}

  SimulatedObjectiveCVaR(const ROL::Ptr<SampleGenerator<Real> > & sampler,
                         const ROL::Ptr<Objective_SimOpt<Real> > & pobj,
                         const ROL::Ptr<PlusFunction<Real> > & pfunc,
                         const Real & alpha)
    : sampler_(sampler), pobj_(pobj), pfunc_(pfunc), alpha_(alpha) {}

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    pobj_->update(x,flag,iter);
  }

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) {
    pobj_->update(x,type,iter);
  }

  Real value(const Vector<Real> &x,
             Real &tol) {
    const Vector_SimOpt<Real> &uz = dynamic_cast<const Vector_SimOpt<Real>&>(x);
    ROL::Ptr<const Vector<Real> > uptr = uz.get_1();
    ROL::Ptr<const Vector<Real> > zptr = uz.get_2();
    const SimulatedVector<Real> &pu = dynamic_cast<const SimulatedVector<Real>&>(*uptr);
    const RiskVector<Real> &rz = dynamic_cast<const RiskVector<Real>&>(*zptr);
    Real t = (*rz.getStatistic(0))[0];
    ROL::Ptr<const Vector<Real> > z = rz.getVector();

    std::vector<Real> param;
    Real weight(0), one(1);
    Real val     = 0;
    Real tmpval  = 0;
    Real tmpsum  = 0;
    Real tmpplus = 0;
    for (typename std::vector<SimulatedVector<Real> >::size_type i=0; i<pu.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pobj_->setParameter(param);
      //tmpval = pobj_->value(*(pu.get(i)), *zptr, tol);
      pobj_->update(*(pu.get(i)), *z);
      tmpval = pobj_->value(*(pu.get(i)), *z, tol);
      tmpplus = pfunc_->evaluate(tmpval-t, 0);
      tmpsum += tmpplus*weight;
    }
    sampler_->sumAll(&tmpsum, &val, 1);
    val *= (one/(one-alpha_));
    val += t;
    return val;
  }
 
  virtual void gradient(Vector<Real> &g,
                        const Vector<Real> &x,
                        Real &tol) {
    g.zero();
    // split x
    const Vector_SimOpt<Real> &xuz = dynamic_cast<const Vector_SimOpt<Real>&>(x);
    ROL::Ptr<const Vector<Real> > xuptr = xuz.get_1();
    ROL::Ptr<const Vector<Real> > xzptr = xuz.get_2();
    const SimulatedVector<Real> &pxu = dynamic_cast<const SimulatedVector<Real>&>(*xuptr);
    const RiskVector<Real> &rxz = dynamic_cast<const RiskVector<Real>&>(*xzptr);
    Real xt = (*rxz.getStatistic(0))[0];
    ROL::Ptr<const Vector<Real> > xz = rxz.getVector();
    // split g
    Vector_SimOpt<Real> &guz = dynamic_cast<Vector_SimOpt<Real>&>(g);
    ROL::Ptr<Vector<Real> > guptr = guz.get_1();
    ROL::Ptr<Vector<Real> > gzptr = guz.get_2();
    SimulatedVector<Real> &pgu = dynamic_cast<SimulatedVector<Real>&>(*guptr);
    RiskVector<Real> &rgz = dynamic_cast<RiskVector<Real>&>(*gzptr);
    ROL::Ptr<Vector<Real> > gz = rgz.getVector();

    std::vector<Real> param;
    Real weight(0), one(1), sum(0), tmpsum(0), tmpval(0), tmpplus(0);
    //ROL::Ptr<Vector<Real> > tmp1 = gzptr->clone();
    //ROL::Ptr<Vector<Real> > tmp2 = gzptr->clone();
    ROL::Ptr<Vector<Real> > tmp1 = gz->clone();
    ROL::Ptr<Vector<Real> > tmp2 = gz->clone();
    for (typename std::vector<SimulatedVector<Real> >::size_type i=0; i<pgu.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pobj_->setParameter(param);
      pobj_->update(*(pxu.get(i)), *xz);
      //tmpval = pobj_->value(*(pxu.get(i)), *xzptr, tol);
      tmpval = pobj_->value(*(pxu.get(i)), *xz, tol);
      tmpplus = pfunc_->evaluate(tmpval-xt, 1);
      tmpsum += weight*tmpplus;
      //Vector_SimOpt<Real> xi(ROL::constPtrCast<Vector<Real> >(pxu.get(i)), ROL::constPtrCast<Vector<Real> >(xzptr));
      Vector_SimOpt<Real> xi(ROL::constPtrCast<Vector<Real> >(pxu.get(i)), ROL::constPtrCast<Vector<Real> >(xz));
      Vector_SimOpt<Real> gi(pgu.get(i), tmp1);
      pobj_->gradient(gi, xi, tol);
      gi.scale(weight*tmpplus);
      tmp2->plus(*tmp1);
      pgu.get(i)->scale(one/(one-alpha_));
    }
    //sampler_->sumAll(*tmp2, *gzptr);
    //gzptr->scale(one/(one-alpha_));
    sampler_->sumAll(*tmp2, *gz);
    gz->scale(one/(one-alpha_));
    sampler_->sumAll(&tmpsum, &sum, 1);
    rgz.setStatistic(one - (one/(one-alpha_))*sum,0);
  }

/*
  virtual void hessVec(Vector<Real> &hv,
                       const Vector<Real> &v,
                       const Vector<Real> &x,
                       Real &tol) {
    hv.zero();
    // split x
    const Vector_SimOpt<Real> &xuz = dynamic_cast<const Vector_SimOpt<Real>&>(x);
    ROL::Ptr<const Vector<Real> > xuptr = xuz.get_1();
    ROL::Ptr<const Vector<Real> > xzptr = xuz.get_2();
    const SimulatedVector<Real> &pxu = dynamic_cast<const SimulatedVector<Real>&>(*xuptr);
    // split v
    const Vector_SimOpt<Real> &vuz = dynamic_cast<const Vector_SimOpt<Real>&>(v);
    ROL::Ptr<const Vector<Real> > vuptr = vuz.get_1();
    ROL::Ptr<const Vector<Real> > vzptr = vuz.get_2();
    const SimulatedVector<Real> &pvu = dynamic_cast<const SimulatedVector<Real>&>(*vuptr);
    // split hv
    Vector_SimOpt<Real> &hvuz = dynamic_cast<Vector_SimOpt<Real>&>(hv);
    ROL::Ptr<Vector<Real> > hvuptr = hvuz.get_1();
    ROL::Ptr<Vector<Real> > hvzptr = hvuz.get_2();
    SimulatedVector<Real> &phvu = dynamic_cast<SimulatedVector<Real>&>(*hvuptr);

    std::vector<Real> param;
    Real weight(0);
    ROL::Ptr<Vector<Real> > tmp1 = hvzptr->clone();
    ROL::Ptr<Vector<Real> > tmp2 = hvzptr->clone();
    for (typename std::vector<SimulatedVector<Real> >::size_type i=0; i<phvu.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pobj_->setParameter(param);
      Vector_SimOpt<Real> xi(ROL::constPtrCast<Vector<Real> >(pxu.get(i)), ROL::constPtrCast<Vector<Real> >(xzptr));
      Vector_SimOpt<Real> vi(ROL::constPtrCast<Vector<Real> >(pvu.get(i)), ROL::constPtrCast<Vector<Real> >(vzptr));
      Vector_SimOpt<Real> hvi(phvu.get(i), tmp1);
      pobj_->update(xi);
      pobj_->hessVec(hvi, vi, xi, tol);
      hvi.scale(weight);
      tmp2->plus(*tmp1);
    }
    sampler_->sumAll(*tmp2, *hvzptr);
  }
*/

}; // class SimulatedObjective

} // namespace ROL

#endif
