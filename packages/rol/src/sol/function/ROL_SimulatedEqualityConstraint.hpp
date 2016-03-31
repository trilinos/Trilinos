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
  const Teuchos::RCP<ParametrizedEqualityConstraint_SimOpt> > pcon_;

public:

  virtual ~SimulatedEqualityConstraint() {}

  SimulatedEqualityConstraint(const Teuchos::RCP<SampleGenerator<Real> > & sampler,
                              const Teuchos::RCP<ParametrizedEqualityConstraint_SimOpt> > & pcon)
    : sampler_(sampler), pcon_(pcon) {}

  void value(Vector<Real> &c,
             const Vector<Real> &x,
             Real &tol) {
    SimulatedVector<Real> &pc = Teuchos::dyn_cast<SimulatedVector<Real> >(c);
    const Vector_SimOpt<Real> &uz = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(x);
    Teuchos::RCP<const Vector<Real> > uptr = uz.get_1();
    Teuchos::RCP<const Vector<Real> > zptr = uz.get_2();
    const SimulatedVector<Real> &pu = Teuchos::dyn_cast<const SimulatedVector<Real> >(*uptr);

    std::vector<Real> param;
    Real weight(0);
    for (std::vector<SimulatedVector<Real> >::size_type i=0; i<pu.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pcon_->setParameter(param);
      pcon_->value(*(pc.get(i)), *(pu.get(i)), *zptr, tol);
      pc.get(i)->scale(weight);
    }

  }

 
  virtual void applyJacobian(Vector<Real> &jv,
                             const Vector<Real> &v,
                             const Vector<Real> &x,
                             Real &tol) {
    SimulatedVector<Real> &pjv = Teuchos::dyn_cast<SimulatedVector<Real> >(jv);
    // split x
    const Vector_SimOpt<Real> &xuz = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(x);
    Teuchos::RCP<Vector<Real> > xuptr = xuz.get_1();
    Teuchos::RCP<Vector<Real> > xzptr = xuz.get_2();
    Vector<Real> &xu = *xuptr;
    Vector<Real> &xz = *xzptr;
    SimulatedVector<Real> &pxu = Teuchos::dyn_cast<SimulatedVector<Real> >(xu);
    // split v
    const Vector_SimOpt<Real> &vuz = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(v);
    Teuchos::RCP<Vector<Real> > vuptr = vuz.get_1();
    Teuchos::RCP<Vector<Real> > vzptr = vuz.get_2();
    Vector<Real> &vu = *vuptr;
    Vector<Real> &vz = *vzptr;
    SimulatedVector<Real> &pvu = Teuchos::dyn_cast<SimulatedVector<Real> >(vu);

    std::vector<Real> param;
    Real weight(0);
    for (std::vector<SimulatedVector<Real> >::size_type i=0; i<pvu.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pcon_->setParameter(param);
      Vector_SimOpt<Real> xi(pxu.get(i), xzptr);
      Vector_SimOpt<Real> vi(pvu.get(i), vzptr);
      pcon_->applyJacobian(*(pjv.get(i)), vi, xi, tol);
      pjv.get(i)->scale(weight);
    }
    
  }

}; // class SimulatedEqualityConstraint

} // namespace ROL

#endif
