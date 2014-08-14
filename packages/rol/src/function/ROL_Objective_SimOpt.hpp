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

#ifndef ROL_OBJECTIVE_SIMOPT_H
#define ROL_OBJECTIVE_SIMOPT_H

#include "ROL_Objective.hpp"
#include "ROL_EqualityConstraint_SimOpt.hpp"
#include "ROL_Vector_SimOpt.hpp"

/** \class ROL::Objective_SimOpt
    \brief Provides the interface to evaluate simulation-based reduced objective functions.
*/


namespace ROL {

template <class Real>
class Objective_SimOpt : public Objective<Real> {
private:
  Teuchos::RCP<Objective<Real> > obj_;
  Teuchos::RCP<EqualityConstraint_SimOpt<Real> > con_;
  Teuchos::RCP<Vector<Real> > state_;
  Teuchos::RCP<Vector<Real> > adjoint_;

  bool is_state_computed_;
  bool is_adjoint_computed_;

public:
  Objective_SimOpt(Teuchos::RCP<Objective<Real> > &obj, Teuchos::RCP<EqualityConstraint_SimOpt<Real> > &con, 
                   Teuchos::RCP<Vector<Real> > &state, Teuchos::RCP<Vector<Real> > &adjoint) 
    : obj_(obj), con_(con), state_(state), adjoint_(adjoint) {
    is_state_computed_   = false;
    is_adjoint_computed_ = false;
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    // Solve state equation at new iterate
    Real tol = std::sqrt(ROL_EPSILON);
    Teuchos::RCP<Vector<Real> > z = x.clone(); z->set(x);
    (this->con_)->solve(*(this->state_),*z,tol);
    Vector_SimOpt<Real> xs(this->state_,z);
    // Update full objective function
    (this->obj_)->update(xs,flag,iter);
    // Update equality constraint
    (this->con_)->update(xs,flag,iter);
    // Reset storage flags
    this->is_state_computed_ = true;
    this->is_adjoint_computed_ = false;
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    // Solve state equation if not done already
    Teuchos::RCP<Vector<Real> > z = x.clone(); z->set(x);
    if (!this->is_state_computed_) {
      (this->con_)->solve(*(this->state_),*z,tol);
    }
    // Build a state/control vector
    Vector_SimOpt<Real> xs(this->state_,z);
    // Get objective function value
    return (this->obj_)->value(xs,tol);
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    // Solve state equation if not done already
    Teuchos::RCP<Vector<Real> > z = x.clone(); z->set(x);
    if(!this->is_state_computed_) {
      (this->con_)->solve(*(this->state_),*z,tol);
    }
    // Build a state/control vector
    Vector_SimOpt<Real> xs(this->state_,z);
    // Evaluate the full gradient
    Teuchos::RCP<Vector<Real> > gu = (this->state_)->clone();
    Teuchos::RCP<Vector<Real> > gz = x.clone();
    Vector_SimOpt<Real> gs(gu,gz);
    (this->obj_)->gradient(gs,xs,tol);
    // Solve adjoint equation if not done already
    if(!this->is_adjoint_computed_) {
      (this->con_)->applyInverseAdjointJacobian_1(*(this->adjoint_),*(gs.get_1()),*(this->state_),*z,tol);
      (this->adjoint_)->scale(-1.0);
    }
    // Build gradient
    (this->con_)->applyAdjointJacobian_2(g,*(this->adjoint_),*(this->state_),*z,tol);
    g.plus(*(gs.get_2()));
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    // Solve state equation if not done already
    Teuchos::RCP<Vector<Real> > z = x.clone(); z->set(x);
    if(!this->is_state_computed_) {
      (this->con_)->solve(*(this->state_),*z,tol);
    }
    // Build a state/control vector
    Vector_SimOpt<Real> xs(this->state_,z);
    // Solve adjoint equation if not done already
    if(!this->is_adjoint_computed_) {
      // Evaluate the full gradient
      Teuchos::RCP<Vector<Real> > gu = (this->state_)->clone();
      Teuchos::RCP<Vector<Real> > gz = x.clone();
      Vector_SimOpt<Real> gs(gu,gz);
      (this->obj_)->gradient(gs,xs,tol);
      // Solve adjoint equation
      (this->con_)->applyInverseAdjointJacobian_1(*(this->adjoint_),*(gs.get_1()),*(this->state_),*z,tol);
      (this->adjoint_)->scale(-1.0);
    }
    // Solve state sensitivity equation
    Teuchos::RCP<Vector<Real> > Bv = (this->state_)->clone();
    (this->con_)->applyJacobian_2(*Bv,v,*(this->state_),*z,tol);
    Bv->scale(-1.0);
    Teuchos::RCP<Vector<Real> > s = (this->state_)->clone();
    (this->con_)->applyInverseJacobian_1(*s,*Bv,*(this->state_),*z,tol);
    // Evaluate full hessVec in the direction (s,v)
    Teuchos::RCP<Vector<Real> > vp = v.clone(); vp->set(v);
    Vector_SimOpt<Real> vs(s,vp);
    Teuchos::RCP<Vector<Real> > hvu = (this->state_)->clone();
    Teuchos::RCP<Vector<Real> > hvz = x.clone();
    Vector_SimOpt<Real> hvs(hvu,hvz);
    (this->obj_)->hessVec(hvs,vs,xs,tol);
    // Apply adjoint Hessian of constraint
    Teuchos::RCP<Vector<Real> > hcu = (this->state_)->clone();
    Teuchos::RCP<Vector<Real> > hcz = x.clone();
    Vector_SimOpt<Real> hcs(hcu,hcz);
    (this->con_)->applyAdjointHessian(hcs,*(this->adjoint_),vs,xs,tol);
    // Solve adjoint sensitivity equation
    Teuchos::RCP<Vector<Real> > r = (this->state_)->clone();
    r->set(*(hvs.get_1()));
    r->plus(*(hcs.get_1())); 
    r->scale(-1.0);
    Teuchos::RCP<Vector<Real> > p = (this->adjoint_)->clone();
    (this->con_)->applyInverseAdjointJacobian_1(*p,*r,*(this->state_),*z,tol);
    // Build hessVec
    (this->con_)->applyAdjointJacobian_2(hv,*p,*(this->state_),*z,tol);
    hv.plus(*(hcs.get_2()));
    hv.plus(*(hvs.get_2()));
  }

  virtual void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    Pv.set(v);
  }

}; // class Step

} // namespace ROL

#endif
