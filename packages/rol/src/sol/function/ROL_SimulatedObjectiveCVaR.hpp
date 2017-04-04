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
  const Teuchos::RCP<SampleGenerator<Real> > sampler_;
  const Teuchos::RCP<Objective_SimOpt<Real> > pobj_;
  const Teuchos::RCP<PlusFunction<Real> > pfunc_;
  const Real alpha_;

public:

  virtual ~SimulatedObjectiveCVaR() {}

  SimulatedObjectiveCVaR(const Teuchos::RCP<SampleGenerator<Real> > & sampler,
                         const Teuchos::RCP<Objective_SimOpt<Real> > & pobj,
                         const Teuchos::RCP<PlusFunction<Real> > & pfunc,
                         const Real & alpha)
    : sampler_(sampler), pobj_(pobj), pfunc_(pfunc), alpha_(alpha) {}

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    pobj_->update(x,flag,iter);
  }

  Real value(const Vector<Real> &x,
             Real &tol) {
    const Vector_SimOpt<Real> &uz = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(x);
    Teuchos::RCP<const Vector<Real> > uptr = uz.get_1();
    Teuchos::RCP<const Vector<Real> > zptr = uz.get_2();
    const SimulatedVector<Real> &pu = Teuchos::dyn_cast<const SimulatedVector<Real> >(*uptr);
    const RiskVector<Real> &rz = Teuchos::dyn_cast<const RiskVector<Real> >(*zptr);
    Real t = rz.getStatistic(0);
    Teuchos::RCP<const Vector<Real> > z = rz.getVector();

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
    const Vector_SimOpt<Real> &xuz = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(x);
    Teuchos::RCP<const Vector<Real> > xuptr = xuz.get_1();
    Teuchos::RCP<const Vector<Real> > xzptr = xuz.get_2();
    const SimulatedVector<Real> &pxu = Teuchos::dyn_cast<const SimulatedVector<Real> >(*xuptr);
    const RiskVector<Real> &rxz = Teuchos::dyn_cast<const RiskVector<Real> >(*xzptr);
    Real xt = rxz.getStatistic(0);
    Teuchos::RCP<const Vector<Real> > xz = rxz.getVector();
    // split g
    Vector_SimOpt<Real> &guz = Teuchos::dyn_cast<Vector_SimOpt<Real> >(g);
    Teuchos::RCP<Vector<Real> > guptr = guz.get_1();
    Teuchos::RCP<Vector<Real> > gzptr = guz.get_2();
    SimulatedVector<Real> &pgu = Teuchos::dyn_cast<SimulatedVector<Real> >(*guptr);
    RiskVector<Real> &rgz = Teuchos::dyn_cast<RiskVector<Real> >(*gzptr);
    Teuchos::RCP<Vector<Real> > gz = rgz.getVector();

    std::vector<Real> param;
    Real weight(0), one(1), sum(0), tmpsum(0), tmpval(0), tmpplus(0);
    //Teuchos::RCP<Vector<Real> > tmp1 = gzptr->clone();
    //Teuchos::RCP<Vector<Real> > tmp2 = gzptr->clone();
    Teuchos::RCP<Vector<Real> > tmp1 = gz->clone();
    Teuchos::RCP<Vector<Real> > tmp2 = gz->clone();
    for (typename std::vector<SimulatedVector<Real> >::size_type i=0; i<pgu.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pobj_->setParameter(param);
      pobj_->update(*(pxu.get(i)), *xz);
      //tmpval = pobj_->value(*(pxu.get(i)), *xzptr, tol);
      tmpval = pobj_->value(*(pxu.get(i)), *xz, tol);
      tmpplus = pfunc_->evaluate(tmpval-xt, 1);
      tmpsum += weight*tmpplus;
      //Vector_SimOpt<Real> xi(Teuchos::rcp_const_cast<Vector<Real> >(pxu.get(i)), Teuchos::rcp_const_cast<Vector<Real> >(xzptr));
      Vector_SimOpt<Real> xi(Teuchos::rcp_const_cast<Vector<Real> >(pxu.get(i)), Teuchos::rcp_const_cast<Vector<Real> >(xz));
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
    rgz.setStatistic(one - (one/(one-alpha_))*sum);
  }

/*
  virtual void hessVec(Vector<Real> &hv,
                       const Vector<Real> &v,
                       const Vector<Real> &x,
                       Real &tol) {
    hv.zero();
    // split x
    const Vector_SimOpt<Real> &xuz = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(x);
    Teuchos::RCP<const Vector<Real> > xuptr = xuz.get_1();
    Teuchos::RCP<const Vector<Real> > xzptr = xuz.get_2();
    const SimulatedVector<Real> &pxu = Teuchos::dyn_cast<const SimulatedVector<Real> >(*xuptr);
    // split v
    const Vector_SimOpt<Real> &vuz = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(v);
    Teuchos::RCP<const Vector<Real> > vuptr = vuz.get_1();
    Teuchos::RCP<const Vector<Real> > vzptr = vuz.get_2();
    const SimulatedVector<Real> &pvu = Teuchos::dyn_cast<const SimulatedVector<Real> >(*vuptr);
    // split hv
    Vector_SimOpt<Real> &hvuz = Teuchos::dyn_cast<Vector_SimOpt<Real> >(hv);
    Teuchos::RCP<Vector<Real> > hvuptr = hvuz.get_1();
    Teuchos::RCP<Vector<Real> > hvzptr = hvuz.get_2();
    SimulatedVector<Real> &phvu = Teuchos::dyn_cast<SimulatedVector<Real> >(*hvuptr);

    std::vector<Real> param;
    Real weight(0);
    Teuchos::RCP<Vector<Real> > tmp1 = hvzptr->clone();
    Teuchos::RCP<Vector<Real> > tmp2 = hvzptr->clone();
    for (typename std::vector<SimulatedVector<Real> >::size_type i=0; i<phvu.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pobj_->setParameter(param);
      Vector_SimOpt<Real> xi(Teuchos::rcp_const_cast<Vector<Real> >(pxu.get(i)), Teuchos::rcp_const_cast<Vector<Real> >(xzptr));
      Vector_SimOpt<Real> vi(Teuchos::rcp_const_cast<Vector<Real> >(pvu.get(i)), Teuchos::rcp_const_cast<Vector<Real> >(vzptr));
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
