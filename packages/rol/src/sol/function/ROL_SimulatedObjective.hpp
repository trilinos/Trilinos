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

#ifndef ROL_SIMULATED_OBJECTIVE_H
#define ROL_SIMULATED_OBJECTIVE_H

#include "ROL_SimulatedVector.hpp"
#include "ROL_Objective_SimOpt.hpp"

namespace ROL {

template <class Real>
class SimulatedObjective : public Objective<Real> {
private:
  const ROL::Ptr<SampleGenerator<Real> > sampler_;
  const ROL::Ptr<Objective_SimOpt<Real> > pobj_;

public:

  virtual ~SimulatedObjective() {}

  SimulatedObjective(const ROL::Ptr<SampleGenerator<Real> > & sampler,
                     const ROL::Ptr<Objective_SimOpt<Real> > & pobj)
    : sampler_(sampler), pobj_(pobj) {}

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {}

  Real value(const Vector<Real> &x,
             Real &tol) {
    const Vector_SimOpt<Real> &uz = dynamic_cast<const Vector_SimOpt<Real>&>(x);
    ROL::Ptr<const Vector<Real> > uptr = uz.get_1();
    ROL::Ptr<const Vector<Real> > zptr = uz.get_2();
    const SimulatedVector<Real> &pu = dynamic_cast<const SimulatedVector<Real>&>(*uptr);

    std::vector<Real> param;
    Real weight(0);
    Real val = 0;
    Real tmpval = 0;
    Real tmpsum = 0;
    for (typename std::vector<SimulatedVector<Real> >::size_type i=0; i<pu.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pobj_->setParameter(param);
      pobj_->update(*(pu.get(i)), *zptr);
      tmpval = pobj_->value(*(pu.get(i)), *zptr, tol);
      tmpsum += tmpval*weight;
    }
    sampler_->sumAll(&tmpsum, &val, 1);
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
    // split g
    Vector_SimOpt<Real> &guz = dynamic_cast<Vector_SimOpt<Real>&>(g);
    ROL::Ptr<Vector<Real> > guptr = guz.get_1();
    ROL::Ptr<Vector<Real> > gzptr = guz.get_2();
    SimulatedVector<Real> &pgu = dynamic_cast<SimulatedVector<Real>&>(*guptr);

    std::vector<Real> param;
    Real weight(0);
    ROL::Ptr<Vector<Real> > tmp1 = gzptr->clone();
    ROL::Ptr<Vector<Real> > tmp2 = gzptr->clone();
    for (typename std::vector<SimulatedVector<Real> >::size_type i=0; i<pgu.numVectors(); ++i) {
      param = sampler_->getMyPoint(static_cast<int>(i));
      weight = sampler_->getMyWeight(static_cast<int>(i));
      pobj_->setParameter(param);
      Vector_SimOpt<Real> xi(ROL::constPtrCast<Vector<Real> >(pxu.get(i)), ROL::constPtrCast<Vector<Real> >(xzptr));
      Vector_SimOpt<Real> gi(pgu.get(i), tmp1);
      pobj_->update(xi);
      pobj_->gradient(gi, xi, tol);
      gi.scale(weight);
      tmp2->plus(*tmp1);
    }
    sampler_->sumAll(*tmp2, *gzptr);

  }


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

}; // class SimulatedObjective

} // namespace ROL

#endif
