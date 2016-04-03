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

#ifndef ROL_SIMULATED_EQUALITY_CONSTRAINT_H
#define ROL_SIMULATED_EQUALITY_CONSTRAINT_H

#include "ROL_SimulatedVector.hpp"
#include "ROL_ParametrizedEqualityConstraint_SimOpt.hpp"

namespace ROL {

template <class Real>
class SimulatedEqualityConstraint : public EqualityConstraint<Real> {
private:
  const Teuchos::RCP<SampleGenerator<Real> > sampler_;
  const Teuchos::RCP<ParametrizedEqualityConstraint_SimOpt<Real> > pcon_;
  const bool useWeights_;

public:

  virtual ~SimulatedEqualityConstraint() {}

  SimulatedEqualityConstraint(const Teuchos::RCP<SampleGenerator<Real> > & sampler,
                              const Teuchos::RCP<ParametrizedEqualityConstraint_SimOpt<Real> > & pcon,
                              const bool useWeights = true)
    : sampler_(sampler), pcon_(pcon), useWeights_(useWeights) {}

  void value(Vector<Real> &c,
             const Vector<Real> &x,
             Real &tol) {
    c.zero();
    SimulatedVector<Real> &pc = Teuchos::dyn_cast<SimulatedVector<Real> >(c);
    const Vector_SimOpt<Real> &uz = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(x);
    Teuchos::RCP<const Vector<Real> > uptr = uz.get_1();
    Teuchos::RCP<const Vector<Real> > zptr = uz.get_2();
    const SimulatedVector<Real> &pu = Teuchos::dyn_cast<const SimulatedVector<Real> >(*uptr);

    std::vector<Real> param;
    Real weight(0), one(1);
    for (typename std::vector<SimulatedVector<Real> >::size_type i=0; i<pu.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pcon_->setParameter(param);
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
    SimulatedVector<Real> &pjv = Teuchos::dyn_cast<SimulatedVector<Real> >(jv);
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

    std::vector<Real> param;
    Real weight(0), one(1);
    for (typename std::vector<SimulatedVector<Real> >::size_type i=0; i<pvu.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pcon_->setParameter(param);
      Vector_SimOpt<Real> xi(Teuchos::rcp_const_cast<Vector<Real> >(pxu.get(i)), Teuchos::rcp_const_cast<Vector<Real> >(xzptr));
      Vector_SimOpt<Real> vi(Teuchos::rcp_const_cast<Vector<Real> >(pvu.get(i)), Teuchos::rcp_const_cast<Vector<Real> >(vzptr));
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
    Vector_SimOpt<Real> &ajvuz = Teuchos::dyn_cast<Vector_SimOpt<Real> >(ajv);
    Teuchos::RCP<Vector<Real> > ajvuptr = ajvuz.get_1();
    Teuchos::RCP<Vector<Real> > ajvzptr = ajvuz.get_2();
    SimulatedVector<Real> &pajvu = Teuchos::dyn_cast<SimulatedVector<Real> >(*ajvuptr);
    // split x
    const Vector_SimOpt<Real> &xuz = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(x);
    Teuchos::RCP<const Vector<Real> > xuptr = xuz.get_1();
    Teuchos::RCP<const Vector<Real> > xzptr = xuz.get_2();
    const SimulatedVector<Real> &pxu = Teuchos::dyn_cast<const SimulatedVector<Real> >(*xuptr);
    // cast v
    const SimulatedVector<Real> &pv = Teuchos::dyn_cast<const SimulatedVector<Real> >(v);

    std::vector<Real> param;
    Real weight(0), one(1);
    Teuchos::RCP<Vector<Real> > tmp1 = ajvzptr->clone();
    Teuchos::RCP<Vector<Real> > tmp2 = ajvzptr->clone();
    for (typename std::vector<SimulatedVector<Real> >::size_type i=0; i<pv.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pcon_->setParameter(param);
      Vector_SimOpt<Real> xi(Teuchos::rcp_const_cast<Vector<Real> >(pxu.get(i)), Teuchos::rcp_const_cast<Vector<Real> >(xzptr));
      Vector_SimOpt<Real> ajvi(pajvu.get(i), tmp1);
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
    Vector_SimOpt<Real> &ahuvuz = Teuchos::dyn_cast<Vector_SimOpt<Real> >(ahuv);
    Teuchos::RCP<Vector<Real> > ahuvuptr = ahuvuz.get_1();
    Teuchos::RCP<Vector<Real> > ahuvzptr = ahuvuz.get_2();
    SimulatedVector<Real> &pahuvu = Teuchos::dyn_cast<SimulatedVector<Real> >(*ahuvuptr);
    // cast u
    const SimulatedVector<Real> &pu = Teuchos::dyn_cast<const SimulatedVector<Real> >(u);
    // split v
    const Vector_SimOpt<Real> &vuz = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(v);
    Teuchos::RCP<const Vector<Real> > vuptr = vuz.get_1();
    Teuchos::RCP<const Vector<Real> > vzptr = vuz.get_2();
    const SimulatedVector<Real> &pvu = Teuchos::dyn_cast<const SimulatedVector<Real> >(*vuptr);
    // split x
    const Vector_SimOpt<Real> &xuz = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(x);
    Teuchos::RCP<const Vector<Real> > xuptr = xuz.get_1();
    Teuchos::RCP<const Vector<Real> > xzptr = xuz.get_2();
    const SimulatedVector<Real> &pxu = Teuchos::dyn_cast<const SimulatedVector<Real> >(*xuptr);

    std::vector<Real> param;
    Real weight(0), one(1);
    Teuchos::RCP<Vector<Real> > tmp1 = ahuvzptr->clone();
    Teuchos::RCP<Vector<Real> > tmp2 = ahuvzptr->clone();
    for (typename std::vector<SimulatedVector<Real> >::size_type i=0; i<pxu.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pcon_->setParameter(param);
      Vector_SimOpt<Real> xi(Teuchos::rcp_const_cast<Vector<Real> >(pxu.get(i)), Teuchos::rcp_const_cast<Vector<Real> >(xzptr));
      Vector_SimOpt<Real> vi(Teuchos::rcp_const_cast<Vector<Real> >(pvu.get(i)), Teuchos::rcp_const_cast<Vector<Real> >(vzptr));
      Vector_SimOpt<Real> ahuvi(pahuvu.get(i), tmp1);
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
    SimulatedVector<Real> &ppv = Teuchos::dyn_cast<SimulatedVector<Real> >(Pv);
    // split x
    const Vector_SimOpt<Real> &xuz = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(x);
    Teuchos::RCP<const Vector<Real> > xuptr = xuz.get_1();
    Teuchos::RCP<const Vector<Real> > xzptr = xuz.get_2();
    const SimulatedVector<Real> &pxu = Teuchos::dyn_cast<const SimulatedVector<Real> >(*xuptr);
    // split g
    const Vector_SimOpt<Real> &guz = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(g);
    Teuchos::RCP<const Vector<Real> > guptr = guz.get_1();
    Teuchos::RCP<const Vector<Real> > gzptr = guz.get_2();
    const SimulatedVector<Real> &pgu = Teuchos::dyn_cast<const SimulatedVector<Real> >(*guptr);
    // cast v
    const SimulatedVector<Real> &pv = Teuchos::dyn_cast<const SimulatedVector<Real> >(v);

    std::vector<Real> param;
    Real weight(0), one(1);
    for (typename std::vector<SimulatedVector<Real> >::size_type i=0; i<pv.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pcon_->setParameter(param);
      Vector_SimOpt<Real> xi(Teuchos::rcp_const_cast<Vector<Real> >(pxu.get(i)), Teuchos::rcp_const_cast<Vector<Real> >(xzptr));
      Vector_SimOpt<Real> gi(Teuchos::rcp_const_cast<Vector<Real> >(pgu.get(i)), Teuchos::rcp_const_cast<Vector<Real> >(gzptr));
      pcon_->applyPreconditioner(*(ppv.get(i)), *(pv.get(i)), xi, gi, tol);
      weight = (useWeights_) ? weight : one;
      ppv.get(i)->scale(one/(weight*weight));
    }

  }


}; // class SimulatedEqualityConstraint

} // namespace ROL

#endif
