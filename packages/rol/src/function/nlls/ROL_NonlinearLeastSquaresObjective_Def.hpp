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

#ifndef ROL_NONLINEARLEASTSQUARESOBJECTIVE_DEF_H
#define ROL_NONLINEARLEASTSQUARESOBJECTIVE_DEF_H

namespace ROL {

template<typename Real>
NonlinearLeastSquaresObjective<Real>::NonlinearLeastSquaresObjective(const Ptr<Constraint<Real>> &con,
                                                                     const Vector<Real> &optvec,
                                                                     const Vector<Real> &convec,
                                                                     const bool GNH)
  : con_(con), GaussNewtonHessian_(GNH) {
  c1_ = convec.clone(); c1dual_ = c1_->dual().clone();
  c2_ = convec.clone();
  x_  = optvec.dual().clone();
}

template<typename Real>
void NonlinearLeastSquaresObjective<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  con_->update(x,type,iter);
  con_->value(*c1_,x,tol);
  c1dual_->set(c1_->dual());
}

template<typename Real>
void NonlinearLeastSquaresObjective<Real>::update( const Vector<Real> &x, bool flag, int iter ) {
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  con_->update(x,flag,iter);
  con_->value(*c1_,x,tol);
  c1dual_->set(c1_->dual());
}

template<typename Real>
Real NonlinearLeastSquaresObjective<Real>::value( const Vector<Real> &x, Real &tol ) {
  Real half(0.5);
  return half*(c1_->dot(*c1_));
}

template<typename Real>
void NonlinearLeastSquaresObjective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  con_->applyAdjointJacobian(g,*c1dual_,x,tol);
}

template<typename Real>
void NonlinearLeastSquaresObjective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  con_->applyJacobian(*c2_,v,x,tol);
  con_->applyAdjointJacobian(hv,c2_->dual(),x,tol);
  if ( !GaussNewtonHessian_ ) {
    con_->applyAdjointHessian(*x_,*c1dual_,v,x,tol);
    hv.plus(*x_);
  }
}

template<typename Real>
void NonlinearLeastSquaresObjective<Real>::precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  con_->applyPreconditioner(Pv,v,x,x.dual(),tol);
}

template<typename Real>
void NonlinearLeastSquaresObjective<Real>::setParameter(const std::vector<Real> &param) {
  Objective<Real>::setParameter(param);
  con_->setParameter(param);
}

} // namespace ROL

#endif
