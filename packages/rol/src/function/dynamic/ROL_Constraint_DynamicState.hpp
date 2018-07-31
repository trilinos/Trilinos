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

#ifndef ROL_CONSTRAINT_DYNAMICSTATE_H
#define ROL_CONSTRAINT_DYNAMICSTATE_H

#include "ROL_DynamicConstraint.hpp"
#include "ROL_Constraint.hpp"

namespace ROL {

template <class Real>
class Constraint_DynamicState : public Constraint<Real> {
private:
  const Ptr<DynamicConstraint<Real>> con_;
  const Ptr<const Vector<Real>>       uo_;
  const Ptr<const Vector<Real>>        z_;
  const Ptr<const TimeStamp<Real>>    ts_;

  Ptr<Vector<Real>> ijv_;
  bool isInit_;

public:
  Constraint_DynamicState(const Ptr<DynamicConstraint<Real>> &con,
                          const Ptr<const Vector<Real>>      &uo,
                          const Ptr<const Vector<Real>>      &z,
                          const Ptr<const TimeStamp<Real>>   &ts)
    : con_(con), uo_(uo), z_(z), ts_(ts), isInit_(false) {}

  void value(Vector<Real> &c,const Vector<Real> &u,Real &tol) {
    con_->value(c,*uo_,u,*z_,*ts_);
  }

  void applyJacobian(Vector<Real> &jv,const Vector<Real> &v,const Vector<Real> &u,Real &tol) {
    con_->applyJacobian_un(jv,v,*uo_,u,*z_,*ts_);
  }

  void applyAdjointJacobian(Vector<Real> &ajv,const Vector<Real> &v,const Vector<Real> &u,Real &tol) {
    con_->applyAdjointJacobian_un(ajv,v,*uo_,u,*z_,*ts_);
  }

  void applyAdjointHessian(Vector<Real> &ahwv,const Vector<Real> &w,const Vector<Real> &v,const Vector<Real> &u,Real &tol) {
    con_->applyAdjointHessian_un_un(ahwv,w,v,*uo_,u,*z_,*ts_);
  }

  void update( const Vector<Real> &u, bool flag = true, int iter = -1 ) {
    //con_->update_un(u,*ts_);
    con_->update(*uo_,u,*z_,*ts_);
  }

  void applyPreconditioner(Vector<Real> &pv,const Vector<Real> &v,const Vector<Real> &u,const Vector<Real> &g,Real &tol) {
    if (!isInit_) {
      ijv_ = u.clone();
      isInit_ = true;
    }
    con_->applyInverseJacobian_un(*ijv_,v,*uo_,u,*z_,*ts_);
    con_->applyInverseAdjointJacobian_un(pv,ijv_->dual(),*uo_,u,*z_,*ts_);
  }

  // Definitions for parametrized (stochastic) equality constraints
  //void setParameter(const std::vector<Real> &param) {
  //  con_->setParameter(param);
  //  Constraint<Real>::setParameter(param);
  //}

}; // class Constraint_State

} // namespace ROL

#endif
