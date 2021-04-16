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

#ifndef ROL_CONSTRAINT_STATE_DEF_H
#define ROL_CONSTRAINT_STATE_DEF_H

namespace ROL {

template<typename Real>
SimConstraint<Real>::SimConstraint(const Ptr<Constraint_SimOpt<Real>> &con,
                                   const Ptr<const Vector<Real>> &z,
                                   bool inSolve) : con_(con), z_(z), inSolve_(inSolve), init_(false) {}

template<typename Real>
void SimConstraint<Real>::update( const Vector<Real> &u, bool flag, int iter ) {
  con_->update_1(u,flag,iter);
  //con_->update(u,*z_,flag,iter);
}

template<typename Real>
void SimConstraint<Real>::update( const Vector<Real> &u, UpdateType type, int iter ) {
  if (inSolve_) con_->solve_update(u,*z_,type,iter);
  else          con_->update_1(u,type,iter);
}

template<typename Real>
void SimConstraint<Real>::value(Vector<Real> &c,const Vector<Real> &u,Real &tol) {
  con_->value(c,u,*z_,tol);
}

template<typename Real>
void SimConstraint<Real>::applyJacobian(Vector<Real> &jv,const Vector<Real> &v,const Vector<Real> &u,Real &tol) {
  con_->applyJacobian_1(jv,v,u,*z_,tol);
}

template<typename Real>
void SimConstraint<Real>::applyAdjointJacobian(Vector<Real> &ajv,const Vector<Real> &v,const Vector<Real> &u,Real &tol) {
  con_->applyAdjointJacobian_1(ajv,v,u,*z_,tol);
}

template<typename Real>
void SimConstraint<Real>::applyAdjointHessian(Vector<Real> &ahwv,const Vector<Real> &w,const Vector<Real> &v,const Vector<Real> &u,Real &tol) {
  con_->applyAdjointHessian_11(ahwv,w,v,u,*z_,tol);
}

template<typename Real>
void SimConstraint<Real>::applyPreconditioner(Vector<Real> &pv,const Vector<Real> &v,const Vector<Real> &u,const Vector<Real> &g,Real &tol) {
  if (!init_) {
    ijv_ = u.clone();
    init_ = true;
  }
  con_->applyInverseJacobian_1(*ijv_,v,u,*z_,tol);
  con_->applyInverseAdjointJacobian_1(pv,ijv_->dual(),u,*z_,tol);
}

template<typename Real>
void SimConstraint<Real>::setParameter(const std::vector<Real> &param) {
  con_->setParameter(param);
  Constraint<Real>::setParameter(param);
}

} // namespace ROL

#endif
