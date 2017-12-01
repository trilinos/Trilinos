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

#ifndef ROL_RISKAVERSE_CONSTRAINT_H
#define ROL_RISKAVERSE_CONSTRAINT_H

#include "ROL_RiskAverseObjective.hpp"
#include "ROL_ConstraintFromObjective.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template <class Real>
class RiskAverseConstraint : public Constraint<Real> {
private:
  ROL::SharedPointer<Objective<Real> > robj_;
  ROL::SharedPointer<Constraint<Real> > con_;
  ROL::SharedPointer<SampleGenerator<Real> > sampler_;

public:
  RiskAverseConstraint(const ROL::SharedPointer<Objective<Real> > &obj,
                       const ROL::SharedPointer<SampleGenerator<Real> > &sampler,
                       ROL::ParameterList &parlist,
                       const int index     = 0)
    : sampler_(sampler) {
    robj_ = ROL::makeShared<RiskAverseObjective<Real>>(obj,parlist,sampler,1,index);
    con_  = ROL::makeShared<ConstraintFromObjective<Real>>(robj_);
  }

  RiskAverseConstraint(const ROL::SharedPointer<Constraint<Real> > &con,
                       const ROL::SharedPointer<SampleGenerator<Real> > &sampler,
                       ROL::ParameterList &parlist,
                       const int index     = 0)
    : sampler_(sampler) {
    try {
      ROL::SharedPointer<ConstraintFromObjective<Real>> cfo
        = ROL::dynamicPointerCast<ConstraintFromObjective<Real>>(con);
      robj_ = ROL::makeShared<RiskAverseObjective<Real>>(cfo->getObjective(),
                parlist,sampler,1,index);
      con_  = ROL::makeShared<ConstraintFromObjective<Real>>(robj_);
    }
    catch (std::exception &e) {
      throw Exception::NotImplemented(">>> ROL::RiskAverseConstraint: Input constraint must be a ConstraintFromObjective!");
    }
  }

  void update(const Vector<Real> &x, bool flag = true, int iter = -1) {
    con_->update(x,flag,iter);
  }

  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
    con_->value(c,x,tol);
  }

  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    con_->applyJacobian(jv,v,x,tol);
  }

  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    con_->applyAdjointJacobian(ajv,v,x,tol);
  }

  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    con_->applyAdjointHessian(ahuv,u,v,x,tol);
  }

}; // class RiskAverseConstraint

} // namespace ROL

#endif
