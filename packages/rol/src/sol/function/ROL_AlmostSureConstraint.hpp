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
  const ROL::Ptr<SampleGenerator<Real> > sampler_;
  const ROL::Ptr<Constraint<Real> >      con_;

  ROL::Ptr<Vector<Real> >                scratch1_;
  ROL::Ptr<Vector<Real> >                scratch2_;
  bool                                       isInitialized_;

public:
  virtual ~AlmostSureConstraint() {}

  AlmostSureConstraint(const ROL::Ptr<SampleGenerator<Real> > &sampler,
                       const ROL::Ptr<Constraint<Real> > &con)
    : sampler_(sampler), con_(con), isInitialized_(false) {}

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {}

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
