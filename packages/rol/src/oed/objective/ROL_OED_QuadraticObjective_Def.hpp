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

#ifndef ROL_OED_QUADRATIC_OBJECTIVE_DEF_HPP
#define ROL_OED_QUADRATIC_OBJECTIVE_DEF_HPP

namespace ROL {
namespace OED {

template<typename Real>
void QuadraticObjective<Real>::initialize(const Vector<Real> &x) {
  if (!isInit_) {
    g_ = x.dual().clone();
    isInit_ = true;
  }
  g_->zero();
}

template<typename Real>
void QuadraticObjective<Real>::apply(const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
  if (!isG_) {
    initialize(u);
    M_->applyJacobian_1(*g_,u,u,z,tol);
    isG_ = true;
  }
}

template<typename Real>
QuadraticObjective<Real>::QuadraticObjective(const Ptr<BilinearConstraint<Real>> &M)
  : ProfiledClass<Real,std::string>("OED::QuadraticObjective"),
    M_(M), isInit_(false), isG_(false) {}

template<typename Real>
Ptr<BilinearConstraint<Real>> QuadraticObjective<Real>::getM() const {
  return M_;
}

template<typename Real>
void QuadraticObjective<Real>::update( const Vector<Real> &u, const Vector<Real> &z, UpdateType type, int iter ) {
  isG_ = false;  
  M_->update(u,z,type,iter);
}

template<typename Real>
Real QuadraticObjective<Real>::value( const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
  startTimer("value");
  apply(u,z,tol);
  Real val = u.dot(g_->dual());
  stopTimer("value");
  return static_cast<Real>(0.5)*val;
}

template<typename Real>
void QuadraticObjective<Real>::gradient_1( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
  startTimer("gradient_1");
  apply(u,z,tol);
  g.set(*g_);
  stopTimer("gradient_1");
}

template<typename Real>
void QuadraticObjective<Real>::gradient_2( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
  startTimer("gradient_2");
  M_->applyAdjointJacobian_2(g,u,u,z,tol);
  g.scale(static_cast<Real>(0.5));
  stopTimer("gradient_2");
}

template<typename Real>
void QuadraticObjective<Real>::hessVec_11( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
  startTimer("hessVec_11");
  M_->applyJacobian_1(hv,v,u,z,tol);
  stopTimer("hessVec_11");
}

template<typename Real>
void QuadraticObjective<Real>::hessVec_12( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
  startTimer("hessVec_12");
  M_->applyAdjointHessian_21(hv,u,v,u,z,tol);
  stopTimer("hessVec_12");
}

template<typename Real>
void QuadraticObjective<Real>::hessVec_21( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
  startTimer("hessVec_21");
  M_->applyAdjointHessian_12(hv,u,v,u,z,tol);
  stopTimer("hessVec_21");
}

template<typename Real>
void QuadraticObjective<Real>::hessVec_22( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
  startTimer("hessVec_22");
  hv.zero();
  stopTimer("hessVec_22");
}

template<typename Real>
void QuadraticObjective<Real>::setParameter( const std::vector<Real> &param ) {
  Objective_SimOpt<Real>::setParameter(param);
  M_->setParameter(param);
}

} // End OED Namespace
} // End ROL Namespace

#endif
