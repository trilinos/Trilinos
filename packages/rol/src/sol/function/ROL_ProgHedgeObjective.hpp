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
//
#ifndef PROGHEDGEOBJECTIVE_H
#define PROGHEDGEOBJECTIVE_H

#include "ROL_Objective.hpp"

/** @ingroup func_group
    \class ROL::ProgHedgeObjective
    \brief Provides the interface for the progressive hedging objective.

    ---
*/
namespace ROL {

template <class Real>
class ProgHedgeObjective : public Objective<Real> {
private:
  const Ptr<Objective<Real>> obj_;
  Ptr<Vector<Real>> xbar_, x0_, w_;
  Real penaltyParam_;

public:

  ProgHedgeObjective(const Ptr<Objective<Real>> &obj,
                     const Vector<Real> &x,
                     const Real penaltyParam) 
    : obj_(obj), penaltyParam_(penaltyParam) {
    xbar_ = x.clone();
    x0_   = x.clone();
    w_    = x.dual().clone();
  }

  void setData(const Vector<Real> &xbar, const Vector<Real> &w, const Real penaltyParam) {
    xbar_->set(xbar);
    w_->set(w);
    penaltyParam_ = penaltyParam;
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    obj_->update(x,flag,iter);
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    const Real half(0.5), one(1);
    Real val = obj_->value(x,tol);
    Real wx  = x.dot(w_->dual());
    x0_->set(x); x0_->axpy(-one,*xbar_);
    Real xx  = x0_->dot(*x0_);
    return val + wx + half*penaltyParam_*xx;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    const Real one(1);
    obj_->gradient(g,x,tol);
    g.plus(*w_);
    x0_->set(x); x0_->axpy(-one,*xbar_);
    g.axpy(penaltyParam_,*x0_);
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    obj_->hessVec(hv,v,x,tol);
    hv.axpy(penaltyParam_,v);
  }

  void setParameter(const std::vector<Real> &param) {
    obj_->setParameter(param);
    Objective<Real>::setParameter(param);
  }

};

}
#endif
