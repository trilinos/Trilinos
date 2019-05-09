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
  const ROL::Ptr<Constraint<Real> >      con_;
  const ROL::Ptr<SampleGenerator<Real> > xsampler_;
  const ROL::Ptr<BatchManager<Real> >    cbman_;

  ROL::Ptr<Vector<Real> > conVec_;
  ROL::Ptr<Vector<Real> > optVec_;

  bool initialized_;

  void init(const Vector<Real> &c, const Vector<Real> &x) {
    if (!initialized_) {
      conVec_ = c.clone();
      optVec_ = x.dual().clone();
      initialized_ = true;
    }
  }

public:
  RiskNeutralConstraint( const ROL::Ptr<Constraint<Real> >      &con,
                         const ROL::Ptr<SampleGenerator<Real> > &xsampler,
                         const ROL::Ptr<BatchManager<Real> >    &cbman)
    : con_(con), xsampler_(xsampler), cbman_(cbman), initialized_(false) {}

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    con_->update(x,flag,iter);
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
